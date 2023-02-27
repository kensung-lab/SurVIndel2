#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <sstream>
#include <unordered_map>

#include "utils.h"
#include "hsr.h"
#include "remapping.h"
#include "consensus.h"
#include "libs/cptl_stl.h"
#include "libs/ssw_cpp.h"
#include "libs/ssw.h"
#include "libs/IntervalTree.h"
#include "sam_utils.h"
#include "stat_tests.h"

config_t config;
stats_t stats;
std::string workdir, complete_bam_fname, reference_fname;
std::mutex del_out_mtx, dup_out_mtx, out_mtx, log_mtx;

chr_seqs_map_t chr_seqs;
contig_map_t contig_map;

std::ofstream flog;

std::mutex bam_pool_mtx;
std::queue<open_samFile_t*> bam_pool;
open_samFile_t* get_bam_reader(std::string bam_fname) {
    open_samFile_t* o;
    bam_pool_mtx.lock();
    if (!bam_pool.empty()) {
        o = bam_pool.front();
        bam_pool.pop();
    } else {
        o = open_samFile(bam_fname.c_str());
        hts_set_fai_filename(o->file, fai_path(reference_fname.c_str()));
    }
    bam_pool_mtx.unlock();
    return o;
}
void release_bam_reader(open_samFile_t* reader) {
    bam_pool_mtx.lock();
    bam_pool.push(reader);
    bam_pool_mtx.unlock();
}

std::mutex mtx;
std::unordered_map<std::string, std::vector<deletion_t*> > deletions_by_chr;
std::unordered_map<std::string, std::vector<duplication_t*> > duplications_by_chr;
std::vector<double> max_allowed_frac_normalized;


struct pair_w_score_t {
    int rc_idx, lc_idx;
    int score;
    double mm_rate;
    bool trimmed;

    pair_w_score_t(int rc_idx, int lc_idx, int score, double mm_rate, bool trimmed) :
    	rc_idx(rc_idx), lc_idx(lc_idx), score(score), mm_rate(mm_rate), trimmed(trimmed) {}
};

indel_t* realign_consensuses(std::string contig_name, consensus_t* rc_consensus, consensus_t* lc_consensus, bool trimmed_consensuses, StripedSmithWaterman::Aligner& aligner) {

	std::string rc_consensus_seq = rc_consensus->consensus, lc_consensus_seq = lc_consensus->consensus;
	if (trimmed_consensuses) {
		rc_consensus_seq = rc_consensus_seq.substr(0, rc_consensus_seq.length()-rc_consensus->lowq_clip_portion);
		lc_consensus_seq = lc_consensus_seq.substr(lc_consensus->lowq_clip_portion);
	}
	suffix_prefix_aln_t sp_aln = aln_suffix_prefix(rc_consensus_seq, lc_consensus_seq, 1, -4);
	std::string joined_consensus = rc_consensus_seq + lc_consensus_seq.substr(sp_aln.overlap);

	char* chr_seq = chr_seqs.get_seq(contig_name);
	hts_pos_t chr_len = chr_seqs.get_len(contig_name);

	hts_pos_t ref_lh_start = rc_consensus->start - config.clip_penalty - 10; // 10 bp tolerance
	hts_pos_t ref_lh_end = rc_consensus->breakpoint + joined_consensus.length();
	if (ref_lh_start < 0) ref_lh_start = 0;
	if (ref_lh_end >= chr_len) ref_lh_end = chr_len-1;
	hts_pos_t ref_lh_len = ref_lh_end - ref_lh_start;

	hts_pos_t ref_rh_start = lc_consensus->breakpoint - joined_consensus.length();
	hts_pos_t ref_rh_end = lc_consensus->end + config.clip_penalty + 10; // 10 bp tolerance
	// lc_consensus->breakpoint + joined_consensus.length()*2;
	if (ref_rh_start < 0) ref_rh_start = 0;
	if (ref_rh_end >= chr_len) ref_rh_end = chr_len-1;
	hts_pos_t ref_rh_len = ref_rh_end - ref_rh_start;

	std::string source;
	if (!rc_consensus->is_hsr && !lc_consensus->is_hsr) {
		source = "2SR";
	} else if (rc_consensus->is_hsr && lc_consensus->is_hsr) {
		source = "2HSR";
	} else if (!rc_consensus->is_hsr && lc_consensus->is_hsr) {
		source = "SR-HSR";
	} else if (rc_consensus->is_hsr && !lc_consensus->is_hsr) {
		source = "HSR-SR";
	}
	indel_t* indel = remap_consensus(joined_consensus, chr_seq, ref_lh_start, ref_lh_len, ref_rh_start, ref_rh_len, aligner,
			lc_consensus, rc_consensus, source);
	if (indel == NULL) return NULL;

	char ref_with_del[100000];
	if (indel->indel_type() == "DEL") {
		strncpy(ref_with_del, chr_seqs.get_seq(contig_name)+ref_lh_start, indel->start-ref_lh_start);
		strncpy(ref_with_del+indel->start-ref_lh_start, chr_seqs.get_seq(contig_name)+indel->end, ref_rh_end-indel->end);
		int ref_with_del_len = (indel->start-ref_lh_start) + (ref_rh_end-indel->end);
		ref_with_del[ref_with_del_len] = '\0';

		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment aln_to_ref;
		std::string padded_joined_consensus = config.clip_penalty_padding() + joined_consensus + config.clip_penalty_padding();
		aligner.Align(padded_joined_consensus.c_str(), ref_with_del, ref_with_del_len, filter, &aln_to_ref, 0);
		if (aln_to_ref.query_begin > 0 || aln_to_ref.query_end < padded_joined_consensus.length()-1) indel->start = indel->end = 0;

		deletion_t* del = (deletion_t*) indel;
		del->remap_boundary_lower = lc_consensus->remap_boundary;
		del->remap_boundary_upper = rc_consensus->remap_boundary;
	} else if (indel->indel_type() == "DUP") {
		duplication_t* dup = (duplication_t*) indel;
		dup->original_start = lc_consensus->breakpoint, dup->original_end = rc_consensus->breakpoint;
	}

	if (indel->len() <= 0) {
		delete indel;
		return NULL;
	}

	return indel;
}

void overlap_read_with_1SR_rc_consensuses(indel_t* indel, bam1_t* read) {
    if (indel->sr_remap_info->strong_supporting_reads > 5) return;

    hts_pos_t remapped_bp = indel->sr_remap_info->remapped_bp;
    if (get_unclipped_start(read) <= remapped_bp &&
        get_unclipped_end(read) >= remapped_bp+config.min_clip_len) {
        indel->sr_remap_info->tried_reads++;

        std::string read_seq = get_sequence(read);
        std::string consensus = indel->rc_consensus->consensus;
        // we need to constrain the overlap - i.e.: we have a read that overlaps with the remapped clip, let's say in
        // positions read[X..X+n]
        // when we do the suffix prefix overlap the read with the consensus, we need to ensure that the clipped part of
        // the consensus is aligning to read[X..X+n]
        int max_overlap = indel->sr_remap_info->remapped_other_end-read->core.pos + get_left_clip_size(read)+1;
        suffix_prefix_aln_t spa = aln_suffix_prefix(consensus, read_seq, 1, -4, config.min_clip_len, max_overlap);

        std::stringstream ss;
//		ss << indel->to_string(contig_name) << " " << bam_get_qname(read) << " " << spa.overlap << " ";
//		ss << spa.mismatches << " " << get_cigar_code(read) << " " << spa.mismatch_rate() << " ";
        if (spa.overlap >= config.min_clip_len) {
			std::pair<int, int> read_mismatches_pair = get_mismatches_prefix(read, spa.overlap);
			int read_mismatches = read_mismatches_pair.first, read_indels = read_mismatches_pair.second;
            int aln_mismatches = spa.mismatches;

            if (spa.mismatch_rate() <= config.max_seq_error) {
            	ss << (is_left_clipped(read, config.min_clip_len) || read_indels > 1 || aln_mismatches <= read_mismatches) << std::endl;
				// Given a read that aligns properly to a clipped cluster consensus,
				// we accept a read as supporting a deletion if either:
				// 1 - the read is clipped
				// 2 - it has at least 2 indels
				// 3 - the alignment to the consensus is not worse than the alignment to the reference
				// For #2, note that we don't allow indels when aligning a read to a consensus - because Illumina is
				// unlikely to add indels as sequencing errors. Having 1 indel compared to the reference is still fine,
				// because we are trying to predict a deletion - that is, by accepting the alignment,
				// we are accepting that there is 1 indel
				if (is_left_clipped(read, config.min_clip_len) || read_indels > 1 || aln_mismatches <= read_mismatches) {
					indel->sr_remap_info->weak_supporting_reads++;

					if (spa.overlap >= indel->rc_consensus->clip_len + config.min_clip_len) {
						indel->sr_remap_info->strong_supporting_reads++;

						if (!spa.mismatches) indel->sr_remap_info->perfect_supporting_reads++;
					}
				}
            }
        }
        log_mtx.lock();
		flog << ss.str() << std::endl;
		log_mtx.unlock();
    }
}

void overlap_read_with_1SR_lc_consensuses(indel_t* indel, bam1_t* read) {
    if (indel->sr_remap_info->strong_supporting_reads > 5) return;

	hts_pos_t remapped_bp = indel->sr_remap_info->remapped_bp;
	if (get_unclipped_start(read) <= remapped_bp-config.min_clip_len &&
        get_unclipped_end(read) >= remapped_bp) {
        indel->sr_remap_info->tried_reads++;

        std::string read_seq = get_sequence(read);
        int max_overlap = bam_endpos(read)-indel->sr_remap_info->remapped_other_end + get_right_clip_size(read);
        suffix_prefix_aln_t spa = aln_suffix_prefix(read_seq, indel->lc_consensus->consensus, 1, -4,
                                                    config.min_clip_len, max_overlap);

        std::stringstream ss;
//		ss << indel->to_string(contig_name) << " " << bam_get_qname(read) << " " << spa.overlap << " ";
//		ss << spa.mismatches << " " << get_cigar_code(read) << " " << spa.mismatch_rate() << " ";
        if (spa.overlap >= config.min_clip_len) {
            std::pair<int, int> read_mismatches_pair = get_mismatches_suffix(read, spa.overlap);
            int read_mismatches = read_mismatches_pair.first, read_indels = read_mismatches_pair.second;

            if (spa.mismatch_rate() <= config.max_seq_error) {
            	ss << (is_right_clipped(read, config.min_clip_len) || read_indels > 1 || spa.mismatches <= read_mismatches) << std::endl;
				if (is_right_clipped(read, config.min_clip_len) || read_indels > 1 || spa.mismatches <= read_mismatches) {
					indel->sr_remap_info->weak_supporting_reads++;

					if (spa.overlap >= indel->lc_consensus->clip_len + config.min_clip_len) {
						indel->sr_remap_info->strong_supporting_reads++;

						if (!spa.mismatches) indel->sr_remap_info->perfect_supporting_reads++;
					}
				}
            }

        }

        log_mtx.lock();
		flog << ss.str() << std::endl;
		log_mtx.unlock();
    }
}


void build_hsr_consensuses(int id, int contig_id, std::string contig_name, std::vector<consensus_t*>& rc_hsr_consensuses,
	std::vector<consensus_t*>& lc_hsr_consensuses) {
    std::string bam_fname = workdir + "/workspace/" + std::to_string(contig_id) + "-HSR.bam";
    open_samFile_t* bam_file = open_samFile(bam_fname, true);

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    mtx.lock();
    std::cout << "Building HSR for " << contig_name << std::endl;
    mtx.unlock();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::deque<bam1_t*> lc_cluster, rc_cluster;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {

        std::pair<int, int> left_and_right_diffs = compute_left_and_right_differences(read, true);

        bool lc_clipped = left_and_right_diffs.first > left_and_right_diffs.second;
        std::deque<bam1_t*>& cluster = lc_clipped ? lc_cluster : rc_cluster;
        std::vector<consensus_t*>& consensuses = lc_clipped ? lc_hsr_consensuses : rc_hsr_consensuses;

        if (cluster.size() >= 3 && bam_endpos(cluster.front())-read->core.pos < read->core.l_qseq/2) {
            // only use HSRs for repeat regions
            // note the tolerance of 50 bp for imprecisions in TRF annotations
			std::vector<bam_redux_t*> cluster_v;
			for (bam1_t* r : cluster) cluster_v.push_back(new bam_redux_t(r));
			consensus_t* consensus = build_full_consensus(contig_id, cluster_v, lc_clipped);
			if (consensus != NULL) {
				consensus->is_hsr = true;
				consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
				consensuses.push_back(consensus);
			}
			for (bam_redux_t* r : cluster_v) delete r;
        }
        while (!cluster.empty() && bam_endpos(cluster.front())-read->core.pos < read->core.l_qseq/2) {
            bam_destroy1(cluster.front());
            cluster.pop_front();
        }
        cluster.push_back(bam_dup1(read));
    }
    if (rc_cluster.size() >= 3) {
        	std::vector<bam_redux_t*> cluster_v;
		for (bam1_t* r : rc_cluster) cluster_v.push_back(new bam_redux_t(r));
		consensus_t* consensus = build_full_consensus(contig_id, cluster_v, false);
		if (consensus != NULL) {
			consensus->is_hsr = true;
			consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
			rc_hsr_consensuses.push_back(consensus);
		}
		for (bam_redux_t* r : cluster_v) delete r;
    }
    if (lc_cluster.size() >= 3) {
		std::vector<bam_redux_t*> cluster_v;
		for (bam1_t* r : lc_cluster) cluster_v.push_back(new bam_redux_t(r));
		consensus_t* consensus = build_full_consensus(contig_id, cluster_v, true);
		if (consensus != NULL) {
			consensus->is_hsr = true;
			consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
			lc_hsr_consensuses.push_back(consensus);
		}
		for (bam_redux_t* r : cluster_v) delete r;
    }

    StripedSmithWaterman::Aligner aligner(2, 2, 4, 1, true);
    StripedSmithWaterman::Filter filter;
    filter_well_aligned_to_ref(chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), rc_hsr_consensuses, aligner, filter);
    filter_well_aligned_to_ref(chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), lc_hsr_consensuses, aligner, filter);

    bam_destroy1(read);

    hts_itr_destroy(iter);
    close_samFile(bam_file);

	select_nonoverlapping_clusters(rc_hsr_consensuses);
    select_nonoverlapping_clusters(lc_hsr_consensuses);
}


void find_indels_from_rc_lc_pairs(std::string contig_name, std::vector<consensus_t*>& rc_consensuses, std::vector<consensus_t*>& lc_consensuses,
		std::vector<deletion_t*>& contig_deletions, std::vector<duplication_t*>& contig_duplications, StripedSmithWaterman::Aligner& aligner) {

	// build interval tree of left-clipped consensuses (for quick search)
	std::vector<Interval<int>> consensus_iv;
	for (int i = 0; i < lc_consensuses.size(); i++) {
		consensus_t* lc_anchor = lc_consensuses[i];
		consensus_iv.push_back(Interval<int>(lc_anchor->breakpoint, lc_anchor->breakpoint+1, i));
	}
	IntervalTree<int> consensus_ivtree = IntervalTree<int>(consensus_iv);

	std::vector<pair_w_score_t> rc_lc_scored_pairs;
	for (int i = 0; i < rc_consensuses.size(); i++) {
		consensus_t* rc_anchor = rc_consensuses[i];

		// let us query the interval tree for left-clipped clusters in compatible positions
		// the compatible positions are calculated based on the mates of the right-clipped cluster
		hts_pos_t query_start, query_end;
		const int MAX_SEARCH_DIST = 10000;
		if (rc_anchor->remap_boundary == INT32_MAX) {
			query_start = rc_anchor->breakpoint - MAX_SEARCH_DIST;
			query_end = rc_anchor->breakpoint + MAX_SEARCH_DIST;
		} else {
			query_start = rc_anchor->remap_boundary - config.max_is;
			query_end = rc_anchor->remap_boundary;
		}
		std::vector<Interval<int>> compatible_lc_idxs = consensus_ivtree.findOverlapping(query_start, query_end);

		std::stringstream log_ss;
		for (auto& iv : compatible_lc_idxs) {
			consensus_t* lc_anchor = lc_consensuses[iv.value];
			if (lc_anchor->max_mapq < config.high_confidence_mapq && rc_anchor->max_mapq < config.high_confidence_mapq) continue;

			int min_overlap = (lc_anchor->is_hsr && rc_anchor->is_hsr) ? 50 : std::min(rc_anchor->clip_len, lc_anchor->clip_len)+config.min_clip_len;
			double max_mm_rate = (lc_anchor->is_hsr && rc_anchor->is_hsr) ? 0 : config.max_seq_error;

			suffix_prefix_aln_t spa = aln_suffix_prefix(rc_anchor->consensus, lc_anchor->consensus, 1, -4, min_overlap);

			log_ss << rc_anchor->name() << " " << rc_anchor->consensus << " " << rc_anchor->breakpoint << std::endl;
			log_ss << lc_anchor->name() << " " << lc_anchor->consensus << " " << lc_anchor->breakpoint << std::endl;
			log_ss << spa.overlap << " " << spa.score << " " << spa.mismatches << std::endl;
			log_ss << spa.mismatch_rate() << " " << max_mm_rate << std::endl;

			if (spa.overlap > 0 && spa.mismatch_rate() <= max_mm_rate) {
				rc_lc_scored_pairs.push_back(pair_w_score_t(i, iv.value, spa.score, spa.mismatch_rate(), false));
				log_ss << "ACCEPTED" << std::endl;
			} else { // trim low quality (i.e., supported by less than 2 reads) bases
				// TODO: investigate if we can use base qualities for this
				std::string rc_consensus_trim = rc_anchor->consensus.substr(0, rc_anchor->consensus.length()-rc_anchor->lowq_clip_portion);
				std::string lc_consensus_trim = lc_anchor->consensus.substr(lc_anchor->lowq_clip_portion);
				log_ss << rc_anchor->name() << " " << rc_consensus_trim << " " << rc_anchor->breakpoint << " TRIMMED" << std::endl;
				log_ss << lc_anchor->name() << " " << lc_consensus_trim << " " << lc_anchor->breakpoint << " TRIMMED" << std::endl;

				suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus_trim, lc_consensus_trim, 1, -4, min_overlap);
				log_ss << spa.overlap << " " << spa.score << " " << spa.mismatches << std::endl;
				log_ss << spa.mismatch_rate() << " " << max_mm_rate << std::endl;
				if (spa.overlap > 0 && spa.mismatch_rate() <= max_mm_rate) {
					rc_lc_scored_pairs.push_back(pair_w_score_t(i, iv.value, spa.score, spa.mismatch_rate(), true));
					log_ss << "ACCEPTED" << std::endl;
				}
			}

			log_mtx.lock();
			flog << log_ss.str() << std::endl;
			log_mtx.unlock();
		}
	}

	std::sort(rc_lc_scored_pairs.begin(), rc_lc_scored_pairs.end(),
	            [](const pair_w_score_t& ps1, const pair_w_score_t& ps2) {return ps1.score > ps2.score;});

	std::vector<bool> used_consensus_rc(rc_consensuses.size(), false), used_consensus_lc(lc_consensuses.size(), false);
	for (pair_w_score_t& ps : rc_lc_scored_pairs) {
		if (used_consensus_rc[ps.rc_idx] || used_consensus_lc[ps.lc_idx]) continue;

		consensus_t* rc_anchor = rc_consensuses[ps.rc_idx];
		consensus_t* lc_anchor = lc_consensuses[ps.lc_idx];

		indel_t* indel = realign_consensuses(contig_name, rc_anchor, lc_anchor, ps.trimmed, aligner);
		if (indel == NULL) continue;
		used_consensus_rc[ps.rc_idx] = used_consensus_lc[ps.lc_idx] = true;

		indel->mm_rate = ps.mm_rate;
		indel->extra_info += rc_anchor->name() + "," + lc_anchor->name() + ",";
		indel->extra_info += std::to_string(rc_anchor->start) + "," + std::to_string(lc_anchor->end) + ",";
		indel->extra_info += std::to_string(indel->rc_anchor_start) + "," + std::to_string(indel->lc_anchor_end);
		if (indel->indel_type() == "DEL") {
			contig_deletions.push_back((deletion_t*) indel);
		} else {
			contig_duplications.push_back((duplication_t*) indel);
		}
	}

	remove_marked_consensuses(rc_consensuses, used_consensus_rc);
	remove_marked_consensuses(lc_consensuses, used_consensus_lc);
}

indel_t* find_indel_from_lc_consensus(consensus_t* consensus, hts_pos_t& remapped_other_end, std::string contig_name,
		StripedSmithWaterman::Aligner& aligner) {
	if (!consensus->left_clipped) return NULL;

	hts_pos_t remap_target_start = consensus->remap_boundary;
	hts_pos_t remap_target_end = consensus->remap_boundary + config.max_is;
	if (remap_target_start == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
		remap_target_start = consensus->breakpoint - config.max_is - 2*consensus->consensus.length();
		remap_target_end = consensus->breakpoint + config.max_is + 2*consensus->consensus.length();
	}
	remap_target_start = std::max(remap_target_start, hts_pos_t(0));
	remap_target_end = std::min(remap_target_end, chr_seqs.get_len(contig_name));

	// do not attempt if reference region has Ns - this is because of in our aligner, Ns will always match
	if (remap_target_start >= remap_target_end ||
		has_Ns(chr_seqs.get_seq(contig_name), remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

	std::string clip = consensus->consensus.substr(0, consensus->clip_len);
	std::string padded_clip = config.clip_penalty_padding() + consensus->consensus.substr(0, consensus->clip_len) + config.clip_penalty_padding();
	if (clip.size() < config.min_clip_len) return NULL;

	StripedSmithWaterman::Alignment aln;
	StripedSmithWaterman::Filter filter;
	char* region_start = chr_seqs.get_seq(contig_name) + remap_target_start;
	aligner.Align(padded_clip.c_str(), region_start, remap_target_end - remap_target_start, filter, &aln, 0);
	int best_aln_score = aln.sw_score;
	int new_start = 0;
	indel_t* smallest_indel = NULL;

	log_mtx.lock();
	flog << "LC REMAP: " << consensus->name() << " " << aln.cigar_string << " " << contig_name << ":" << remap_target_start << "-" << remap_target_end << " ";
	flog << padded_clip.c_str() << " " << (int) consensus->max_mapq << std::endl;
	log_mtx.unlock();

	if (is_left_clipped(aln)) {
		clip = clip.substr(consensus->lowq_clip_portion);
		padded_clip = config.clip_penalty_padding() + clip + config.clip_penalty_padding();
//		aligner.Align(padded_clip.c_str(), region_start, remap_target_end - remap_target_start, filter, &aln, 0);

		log_mtx.lock();
		flog << consensus->name() << " " << aln.cigar_string << " " << contig_name << ":" << remap_target_start << "-" << remap_target_end << " ";
		flog << padded_clip.c_str() << " TRIM" << std::endl;
		log_mtx.unlock();
	}

	hts_pos_t ref_rh_end = consensus->end + EXTRA_SEQ;
	if (ref_rh_end >= chr_seqs.get_len(contig_name)) ref_rh_end = chr_seqs.get_len(contig_name)-1;
	hts_pos_t ref_rh_start = consensus->start - config.max_is;
	if (ref_rh_start < 0) ref_rh_start = 0;

	hts_pos_t temp;
	indel_t* indel = remap_consensus(consensus->consensus, chr_seqs.get_seq(contig_name), remap_target_start, remap_target_end-remap_target_start,
			ref_rh_start, ref_rh_end-ref_rh_start, aligner, consensus, NULL, "1SR_LC", remapped_other_end, temp, false);
	if (indel == NULL) return NULL;
	if (indel->indel_type() == "DUP") {
		duplication_t* dup = (duplication_t*) indel;
		dup->original_start = dup->start, dup->original_end = dup->end;
	}
	indel->extra_info += consensus->name();
	return indel;
}

indel_t* find_indel_from_rc_consensus(consensus_t* consensus, hts_pos_t& remapped_other_end, std::string contig_name,
		StripedSmithWaterman::Aligner& aligner) {
	if (consensus->left_clipped) return NULL;

	hts_pos_t remap_target_end = consensus->remap_boundary;
	hts_pos_t remap_target_start = consensus->remap_boundary - config.max_is;
	if (remap_target_end == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
		remap_target_start = consensus->breakpoint - config.max_is - 2*consensus->consensus.length();
		remap_target_end = consensus->breakpoint + config.max_is + 2*consensus->consensus.length();
	}
	remap_target_start = std::max(remap_target_start, hts_pos_t(0));
	remap_target_end = std::min(remap_target_end, chr_seqs.get_len(contig_name));

	if (remap_target_start >= remap_target_end ||
		has_Ns(chr_seqs.get_seq(contig_name), remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

	std::string clip = consensus->consensus.substr(consensus->anchor_len());
	std::string padded_clip = clip + config.clip_penalty_padding();
	if (clip.size() < config.min_clip_len) return NULL;

	StripedSmithWaterman::Alignment aln;
	StripedSmithWaterman::Filter filter;
	char* region_start = chr_seqs.get_seq(contig_name) + remap_target_start;
	aligner.Align(padded_clip.c_str(), region_start, remap_target_end - remap_target_start, filter, &aln, 0);
	int best_aln_score = aln.sw_score;
	int new_start = 0;
	indel_t* smallest_indel = NULL;

	log_mtx.lock();
	flog << "RC REMAP: " <<  consensus->name() << " " << aln.cigar_string << " " << contig_name << ":" << remap_target_start << "-" << remap_target_end << " ";
	flog << padded_clip.c_str() << " " << (int) consensus->max_mapq << std::endl;
	log_mtx.unlock();

	if (is_right_clipped(aln)) {
		clip = clip.substr(0, clip.length()-consensus->lowq_clip_portion);
		padded_clip = clip + config.clip_penalty_padding();
//		aligner.Align(padded_clip.c_str(), region_start, remap_target_end - remap_target_start, filter, &aln, 0);

		log_mtx.lock();
		flog << consensus->name() << " " << aln.cigar_string << " " << contig_name << ":" << remap_target_start << "-" << remap_target_end << " ";
		flog << padded_clip.c_str() << " TRIM" << std::endl;
		log_mtx.unlock();
	}

	hts_pos_t ref_lh_start = consensus->start - EXTRA_SEQ;
	if (ref_lh_start < 0) ref_lh_start = 0;
	hts_pos_t ref_lh_end = consensus->start + consensus->consensus.length();
	if (ref_lh_end >= chr_seqs.get_len(contig_name)) ref_lh_end = chr_seqs.get_len(contig_name)-1;

	hts_pos_t temp;
	indel_t* indel = remap_consensus(consensus->consensus, chr_seqs.get_seq(contig_name), ref_lh_start, ref_lh_end-ref_lh_start, remap_target_start,
			remap_target_end-remap_target_start, aligner, NULL, consensus, "1SR_RC", temp, remapped_other_end, true);
	if (indel == NULL) return NULL;
	if (indel->indel_type() == "DUP") {
		duplication_t* dup = (duplication_t*) indel;
		dup->original_start = dup->start, dup->original_end = dup->end;
	}
	indel->extra_info += consensus->name();
	return indel;
}

void build_clip_consensuses(int id, int contig_id, std::string contig_name, std::vector<deletion_t*>& deletions,
                            std::vector<duplication_t*>& duplications) {

	std::vector<deletion_t*> contig_deletions;
	std::vector<duplication_t*> contig_duplications;

	StripedSmithWaterman::Aligner aligner(1,4,6,1,true);
	StripedSmithWaterman::Filter filter;

    std::string clip_fname = workdir + "/workspace/" + std::to_string(contig_id) + "-CLIP.bam";
    open_samFile_t* bam_file = open_samFile(clip_fname, true);

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::vector<bam_redux_t*> lc_reads, rc_reads;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_left_clipped(read, config.min_clip_len)) {
            lc_reads.push_back(new bam_redux_t(read));
        }
        if (is_right_clipped(read, config.min_clip_len)) {
            rc_reads.push_back(new bam_redux_t(read));
        }
    }

    bam_destroy1(read);

    hts_itr_destroy(iter);
    close_samFile(bam_file);

    std::vector<consensus_t*> rc_consensuses, lc_consensuses;
    std::vector<bam_redux_t*> curr_candidate_cluster;

    auto lc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->start-r2->start) <= config.max_clipped_pos_dist;};
    if (!lc_reads.empty()) {
        for (bam_redux_t* lc_read : lc_reads) {
            if (!curr_candidate_cluster.empty() &&
                !lc_same_cluster(curr_candidate_cluster[0], lc_read)) { // candidate cluster complete
                consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, true);
                curr_candidate_cluster.clear();
                if (full_consensus != NULL) {
                    lc_consensuses.push_back(full_consensus);
                }
            }
            curr_candidate_cluster.push_back(lc_read);
        }
        // process last cluster
        consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, true);
        curr_candidate_cluster.clear();
        if (full_consensus != NULL) {
        	lc_consensuses.push_back(full_consensus);
        }
    }
    auto rc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->end-r2->end) <= config.max_clipped_pos_dist;};
    sort(rc_reads.begin(), rc_reads.end(), [](bam_redux_t* r1, bam_redux_t* r2) { return r1->end < r2->end; });
    if (!rc_reads.empty()) {
        for (bam_redux_t* rc_read : rc_reads) {
            if (!curr_candidate_cluster.empty() &&
                !rc_same_cluster(curr_candidate_cluster[0], rc_read)) { // candidate cluster complete
                consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, false);
                if (full_consensus != NULL) {
                	rc_consensuses.push_back(full_consensus);
                }
                curr_candidate_cluster.clear();
            }
            curr_candidate_cluster.push_back(rc_read);
        }
        // process last cluster
        consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, false);
        curr_candidate_cluster.clear();
        if (full_consensus != NULL) {
            rc_consensuses.push_back(full_consensus);
        }
    }

    for (bam_redux_t* r : lc_reads) delete r;
    for (bam_redux_t* r : rc_reads) delete r;

    auto consensus_cmp = [](consensus_t* c1, consensus_t* c2) {
		return c1->breakpoint < c2->breakpoint;
	};

    // sort by breakpoint position
    sort(rc_consensuses.begin(), rc_consensuses.end(), consensus_cmp);
    sort(lc_consensuses.begin(), lc_consensuses.end(), consensus_cmp);

    // build HSR consensuses
    std::vector<consensus_t*> rc_hsr_consensuses, lc_hsr_consensuses;
    build_hsr_consensuses(0, contig_id, contig_name, rc_hsr_consensuses, lc_hsr_consensuses);
    remove_hsr_overlapping_clipped(rc_hsr_consensuses, rc_consensuses);
    remove_hsr_overlapping_clipped(lc_hsr_consensuses, lc_consensuses);

    std::vector<consensus_t*> full_consensuses;
	for (consensus_t* c : rc_consensuses) full_consensuses.push_back(c);
	for (consensus_t* c : lc_consensuses) full_consensuses.push_back(c);
	std::sort(full_consensuses.begin(), full_consensuses.end(), consensus_cmp);

    // TODO: remove this logging when not needed anymore
	for (consensus_t* consensus : rc_consensuses) {
		log_mtx.lock();
		flog << consensus->name() << " " << consensus->consensus << " " << consensus->breakpoint << " " << consensus->remap_boundary << std::endl;
		log_mtx.unlock();
	}

	log_mtx.lock();
	flog << "N CONSENSUSES: " << contig_name << " " << rc_consensuses.size() << " " << lc_consensuses.size() << " " << rc_hsr_consensuses.size() << " " << lc_hsr_consensuses.size() << std::endl;
	log_mtx.unlock();

	// find RC-LC pairs
	find_indels_from_rc_lc_pairs(contig_name, rc_consensuses, lc_consensuses, contig_deletions, contig_duplications, aligner);
	find_indels_from_rc_lc_pairs(contig_name, rc_hsr_consensuses, lc_consensuses, contig_deletions, contig_duplications, aligner);
	find_indels_from_rc_lc_pairs(contig_name, rc_consensuses, lc_hsr_consensuses, contig_deletions, contig_duplications, aligner);
	find_indels_from_rc_lc_pairs(contig_name, rc_hsr_consensuses, lc_hsr_consensuses, contig_deletions, contig_duplications, aligner);

    /* == Deal with unpaired consensuses == */

    /* == Find possible indels by remapping the clips of unpaired consensuses == */
    std::vector<indel_t*> unpaired_rc_dels, unpaired_lc_dels;
    std::vector<indel_t*> unpaired_rc_dups, unpaired_lc_dups;
    for (int i = 0; i < rc_consensuses.size(); i++) {
		consensus_t* consensus = rc_consensuses[i];
		hts_pos_t remapped_other_end;
		indel_t* smallest_indel = find_indel_from_rc_consensus(consensus, remapped_other_end, contig_name, aligner);

		if (smallest_indel == NULL) continue;
		if (smallest_indel->indel_type() == "DEL") {
			smallest_indel->sr_remap_info = new sr_remap_info_t(smallest_indel->end, remapped_other_end);
			unpaired_rc_dels.push_back(smallest_indel);
		} else {
			smallest_indel->sr_remap_info = new sr_remap_info_t(smallest_indel->start, remapped_other_end);
			unpaired_rc_dups.push_back(smallest_indel);
		}
	}
    for (int i = 0; i < lc_consensuses.size(); i++) {
        consensus_t* consensus = lc_consensuses[i];
		hts_pos_t remapped_other_end;
		indel_t* smallest_indel = find_indel_from_lc_consensus(consensus, remapped_other_end, contig_name, aligner);

		if (smallest_indel == NULL) continue;
		if (smallest_indel->indel_type() == "DEL") {
			smallest_indel->sr_remap_info = new sr_remap_info_t(smallest_indel->start, remapped_other_end);
			unpaired_lc_dels.push_back(smallest_indel);
		} else {
			smallest_indel->sr_remap_info = new sr_remap_info_t(smallest_indel->end, remapped_other_end);
			unpaired_lc_dups.push_back(smallest_indel);
		}
    }

    for (consensus_t* hsr_consensus : rc_hsr_consensuses) {
    	hts_pos_t remapped_other_end;
		indel_t* indel = remap_rc_cluster(hsr_consensus, remapped_other_end, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name),
				aligner);
		if (indel == NULL) continue;
		if (indel->indel_type() == "DEL") {
			indel->sr_remap_info = new sr_remap_info_t(indel->end, remapped_other_end);
			unpaired_rc_dels.push_back(indel);
		} else if (indel->indel_type() == "DUP") {
			indel->sr_remap_info = new sr_remap_info_t(indel->start, remapped_other_end);
			unpaired_rc_dups.push_back(indel);
		}
	}
	for (consensus_t* hsr_consensus : lc_hsr_consensuses) {
		hts_pos_t remapped_other_end;
		indel_t* indel = remap_lc_cluster(hsr_consensus, remapped_other_end, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name),
				aligner);
		if (indel == NULL) continue;
		if (indel->indel_type() == "DEL") {
			indel->sr_remap_info = new sr_remap_info_t(indel->start, remapped_other_end);
			unpaired_lc_dels.push_back(indel);
		} else if (indel->indel_type() == "DUP") {
			indel->sr_remap_info = new sr_remap_info_t(indel->end, remapped_other_end);
			unpaired_lc_dups.push_back(indel);
		}
	}

	auto remapped_bp_cmp = [] (indel_t* i1, indel_t* i2) { return i1->sr_remap_info->remapped_bp < i2->sr_remap_info->remapped_bp; };
	std::sort(unpaired_lc_dels.begin(), unpaired_lc_dels.end(), remapped_bp_cmp);
	std::sort(unpaired_lc_dups.begin(), unpaired_lc_dups.end(), remapped_bp_cmp);
	std::sort(unpaired_rc_dels.begin(), unpaired_rc_dels.end(), remapped_bp_cmp);
	std::sort(unpaired_rc_dups.begin(), unpaired_rc_dups.end(), remapped_bp_cmp);

    bam_file = get_bam_reader(complete_bam_fname);

    // define regions of interest from the remapped clips/breakpoints
    std::vector<char*> regions;
    auto indel_to_region = [&contig_name](const indel_t* indel) {
        std::stringstream ss;
        ss << contig_name << ":" << std::max(hts_pos_t(0), indel->sr_remap_info->remapped_bp-config.read_len) << "-";
        ss << indel->sr_remap_info->remapped_bp+config.read_len;
        char* region = new char[1000];
        strcpy(region, ss.str().c_str());
        return region;
    };
    for (indel_t* indel : unpaired_rc_dels) {
        regions.push_back(indel_to_region(indel));
    }
    for (indel_t* indel : unpaired_lc_dels) {
        regions.push_back(indel_to_region(indel));
    }
    for (indel_t* indel : unpaired_rc_dups) {
        regions.push_back(indel_to_region(indel));
    }
    for (indel_t* indel : unpaired_lc_dups) {
        regions.push_back(indel_to_region(indel));
    }

    if (!regions.empty()) {

    iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    read = bam_init1();

    /* == Find reads that support the remapped breakpoint == */
    out_mtx.lock();
    std::cout << "Computing agreeing reads for " << contig_name << std::endl;
    out_mtx.unlock();
    int curr_rc_del = 0, curr_lc_del = 0, curr_rc_dup = 0, curr_lc_dup = 0;
    int reads_counter = 0;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        reads_counter++;
        if (reads_counter % 10000000 == 0) {
            std::cout << contig_name << " still going!" << std::endl;
        }

        // right-clipped deletions
        while (curr_rc_del < unpaired_rc_dels.size() && (unpaired_rc_dels[curr_rc_del]->sr_remap_info->strong_supporting_reads > 5 ||
        			read->core.pos-read->core.l_qseq > unpaired_rc_dels[curr_rc_del]->sr_remap_info->remapped_bp)) curr_rc_del++;
        for (int i = curr_rc_del; i < unpaired_rc_dels.size() &&
                                  get_unclipped_end(read) > unpaired_rc_dels[i]->sr_remap_info->remapped_bp; i++) {
            if (is_right_clipped(read, config.min_clip_len)) break;

            overlap_read_with_1SR_rc_consensuses(unpaired_rc_dels[i], read);
        }

        // right-clipped duplications
        while (curr_rc_dup < unpaired_rc_dups.size() && (unpaired_rc_dups[curr_rc_dup]->sr_remap_info->strong_supporting_reads > 5 ||
        			read->core.pos-read->core.l_qseq > unpaired_rc_dups[curr_rc_dup]->sr_remap_info->remapped_bp)) curr_rc_dup++;
        for (int i = curr_rc_dup; i < unpaired_rc_dups.size() &&
                                  get_unclipped_end(read) > unpaired_rc_dups[i]->sr_remap_info->remapped_bp; i++) {
            if (is_right_clipped(read, config.min_clip_len)) break;

            overlap_read_with_1SR_rc_consensuses(unpaired_rc_dups[i], read);
        }

        // left-clipped deletions
        while (curr_lc_del < unpaired_lc_dels.size() && (unpaired_lc_dels[curr_lc_del]->sr_remap_info->strong_supporting_reads > 5 ||
        			read->core.pos-read->core.l_qseq > unpaired_lc_dels[curr_lc_del]->sr_remap_info->remapped_bp)) curr_lc_del++;
        for (int i = curr_lc_del; i < unpaired_lc_dels.size() &&
                                  get_unclipped_end(read) > unpaired_lc_dels[i]->sr_remap_info->remapped_bp; i++) {
            if (is_left_clipped(read, config.min_clip_len)) break;

            overlap_read_with_1SR_lc_consensuses(unpaired_lc_dels[i], read);
        }

        // left-clipped duplications
        while (curr_lc_dup < unpaired_lc_dups.size() && (unpaired_lc_dups[curr_lc_dup]->sr_remap_info->strong_supporting_reads > 5 ||
        			read->core.pos-read->core.l_qseq > unpaired_lc_dups[curr_lc_dup]->sr_remap_info->remapped_bp)) curr_lc_dup++;
        for (int i = curr_lc_dup; i < unpaired_lc_dups.size() &&
                                  get_unclipped_end(read) > unpaired_lc_dups[i]->sr_remap_info->remapped_bp; i++) {
            if (is_left_clipped(read, config.min_clip_len)) break;

            overlap_read_with_1SR_lc_consensuses(unpaired_lc_dups[i], read);
        }
    }
    out_mtx.lock();
    std::cout << "Computed agreeing reads for " << contig_name << std::endl;
    out_mtx.unlock();

    for (indel_t* indel : unpaired_rc_dels) {
        if ((indel->sr_remap_info->strong_supporting_reads > 0 && indel->sr_remap_info->weak_supporting_reads >= 3) || indel->rc_consensus->is_hsr || true) {
        	contig_deletions.push_back((deletion_t*) indel);
        }
    }
    for (indel_t* indel : unpaired_lc_dels) {
        if ((indel->sr_remap_info->strong_supporting_reads > 0 && indel->sr_remap_info->weak_supporting_reads >= 3) || indel->lc_consensus->is_hsr || true) {
        	contig_deletions.push_back((deletion_t*) indel);
        }
    }
    for (indel_t* indel : unpaired_rc_dups) {
        if ((indel->sr_remap_info->strong_supporting_reads > 0 && indel->sr_remap_info->weak_supporting_reads >= 3) || indel->rc_consensus->is_hsr || true) {
            contig_duplications.push_back((duplication_t*) indel);
        }
    }
    for (indel_t* indel : unpaired_lc_dups) {
        if ((indel->sr_remap_info->strong_supporting_reads > 0 && indel->sr_remap_info->weak_supporting_reads >= 3) || indel->lc_consensus->is_hsr || true) {
        	contig_duplications.push_back((duplication_t*) indel);
        }
    }

    }
    release_bam_reader(bam_file);

    auto rm_small_indels = [](indel_t* indel) { return indel->len() < config.min_sv_size; };
    del_out_mtx.lock();
	deletions.insert(deletions.end(), contig_deletions.begin(), contig_deletions.end());
    deletions.erase(std::remove_if(deletions.begin(), deletions.end(), rm_small_indels), deletions.end());
	del_out_mtx.unlock();
	dup_out_mtx.lock();
	duplications.insert(duplications.end(), contig_duplications.begin(), contig_duplications.end());
    duplications.erase(std::remove_if(duplications.begin(), duplications.end(), rm_small_indels), duplications.end());
	dup_out_mtx.unlock();
}


void build_consensuses(int id, int contig_id, std::string contig_name) {
    mtx.lock();
    std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
    std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
    mtx.unlock();
    build_clip_consensuses(id, contig_id, contig_name, deletions, duplications);
}

void size_and_depth_filtering(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    out_mtx.lock();
    std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
    std::cout << "Starting confidence intervals computation for " << contig_name << std::endl;
    out_mtx.unlock();
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::vector<double> temp1;
    std::vector<uint32_t> temp2;
    calculate_confidence_interval_size(contig_name, temp1, temp2, deletions, bam_file, config, stats);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    out_mtx.lock();
    std::cout << "Finished confidence interval computation for " << contig_name << ": " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    std::cout << "Starting del depth computation for " << contig_name << std::endl;
    out_mtx.unlock();
    begin = std::chrono::steady_clock::now();
    depth_filter_del(contig_name, deletions, bam_file, config.min_size_for_depth_filtering, config);
    end = std::chrono::steady_clock::now();
    out_mtx.lock();
    std::cout << "Finished del depth computation for " << contig_name << ": " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    out_mtx.unlock();

    out_mtx.lock();
    std::cout << "Starting dup depth computation for " << contig_name << std::endl;
    std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
    out_mtx.unlock();
    begin = std::chrono::steady_clock::now();
    std::vector<duplication_t*> duplications_w_cleanup, duplications_wo_cleanup;
    for (duplication_t* dup : duplications) {
        if (dup->len() >= config.min_size_for_depth_filtering) duplications_w_cleanup.push_back(dup);
        else duplications_wo_cleanup.push_back(dup);
    }
    depth_filter_dup_w_cleanup(contig_name, duplications_w_cleanup, bam_file, stats, config, max_allowed_frac_normalized, workdir);
    depth_filter_dup(contig_name, duplications_wo_cleanup, bam_file, config.min_size_for_depth_filtering, config);
    end = std::chrono::steady_clock::now();
    out_mtx.lock();
    std::cout << "Finished dup depth computation for " << contig_name << ": " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    out_mtx.unlock();

    release_bam_reader(bam_file);
}


int main(int argc, char* argv[]) {

    complete_bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    reference_fname = argv[3];
    std::string sample_name = argv[4];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt");

    chr_seqs.read_fasta_into_map(reference_fname);

//    flog.open(workdir + "/log");
    flog.open("/dev/null");

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        if (contig_name != "chr7") continue;
//        std::future<void> future = thread_pool.push(build_consensuses, contig_id, contig_name);
//        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

    // create VCF out files
	bcf_hdr_t* out_vcf_header = generate_vcf_header(chr_seqs, sample_name, config, full_cmd_str);
    std::string out_vcf_fname = workdir + "/sr.vcf.gz";
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

    std::ifstream as_dist_fin(workdir + "/as_diff_dist.txt");
    int as_diff; double freq;
    std::unordered_map<int, double> as_diff_frequency_table;
    while (as_dist_fin >> as_diff >> freq) {
        as_diff_frequency_table[as_diff] = freq;
    }
    int max_as_diff = 0;
    for (auto& e : as_diff_frequency_table) {
        if (e.first > max_as_diff) max_as_diff = e.first;
    }
    max_allowed_frac_normalized.resize(max_as_diff+1);
    for (int i = 0; i <= max_as_diff; i++) {
        max_allowed_frac_normalized[i] = as_diff_frequency_table[i]/as_diff_frequency_table[0];
    }

    ctpl::thread_pool thread_pool2(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        if (contig_name != "chr7") continue;
        std::future<void> future = thread_pool2.push(size_and_depth_filtering, contig_name);
        futures.push_back(std::move(future));
    }
    thread_pool2.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

    std::unordered_map<std::string, std::vector<bcf1_t*>> bcf_entries;
    bcf1_t* bcf_entry = bcf_init();
    int del_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!deletions_by_chr.count(contig_name)) continue;
		std::vector<deletion_t*>& dels = deletions_by_chr[contig_name];
		std::sort(dels.begin(), dels.end(), [](deletion_t* del1, deletion_t* del2) {
			return std::tie(del1->start, del1->end) < std::tie(del2->start, del2->end);
		});
        for (deletion_t* del : dels) {
        	del->id = "DEL_SR_" + std::to_string(del_id++);
            std::vector<std::string> filters;

            if (del->len() < config.min_sv_size) {
            	filters.push_back("SMALL");
            }
            if (del->len() < config.max_is) {
            	if (del->len()/2 > del->max_conf_size) filters.push_back("SIZE_FILTER");
            }
            if (del->remap_boundary_lower > del->start) {
                filters.push_back("REMAP_BOUNDARY_FILTER");
            } else if (del->remap_boundary_upper < del->end) {
                filters.push_back("REMAP_BOUNDARY_FILTER");
            }
            if (del->len() >= config.min_size_for_depth_filtering &&
                    (del->med_left_flanking_cov*0.74<=del->med_indel_left_cov || del->med_right_flanking_cov*0.74<=del->med_indel_right_cov)) {
            	filters.push_back("DEPTH_FILTER");
            }
            if (del->med_left_flanking_cov > stats.max_depth || del->med_right_flanking_cov > stats.max_depth ||
            	del->med_left_flanking_cov < stats.min_depth || del->med_right_flanking_cov < stats.min_depth ||
				del->med_left_cluster_cov > stats.max_depth || del->med_right_cluster_cov > stats.max_depth) {
                filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (del->med_indel_left_cov > stats.max_depth || del->med_indel_right_cov > stats.max_depth) {
                filters.push_back("ANOMALOUS_DEL_DEPTH");
            }
            if (del->full_junction_score && del->split_junction_score-del->full_junction_score < config.min_score_diff) {
            	filters.push_back("WEAK_SPLIT_ALIGNMENT");
            }
            if ((del->lc_consensus == NULL || del->lc_consensus->max_mapq < config.high_confidence_mapq) &&
            	(del->rc_consensus == NULL || del->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
            }
            if (del->sr_remap_info && (del->sr_remap_info->weak_supporting_reads < 3 || del->sr_remap_info->strong_supporting_reads < 1)) {
            	filters.push_back("WEAK_SUPPORT");
            }

            if (filters.empty()) {
                filters.push_back("PASS");
            }

            del2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, del, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
        }
    }

    int dup_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	if (!duplications_by_chr.count(contig_name)) continue;
    	std::vector<duplication_t*>& dups = duplications_by_chr[contig_name];
    	std::sort(dups.begin(), dups.end(), [](duplication_t* dup1, duplication_t* dup2) {
    		return std::tie(dup1->start, dup1->end) < std::tie(dup2->start, dup2->end);
    	});
        for (duplication_t* dup : dups) {
        	dup->id = "DUP_SR_" + std::to_string(dup_id++);
            std::vector<std::string> filters;

            if (dup->len() >= config.min_size_for_depth_filtering &&
                    (dup->med_left_flanking_cov*1.26>=dup->med_indel_left_cov || dup->med_indel_right_cov<=dup->med_right_flanking_cov*1.26)) {
                // note: using >= so that a 0 0 0 0 depth will not be accepted
                filters.push_back("DEPTH_FILTER");
            }

            if (dup->med_left_flanking_cov > stats.max_depth || dup->med_right_flanking_cov > stats.max_depth ||
                dup->med_left_flanking_cov < stats.min_depth || dup->med_right_flanking_cov < stats.min_depth) {
                filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            }
            if (dup->full_junction_score && dup->split_junction_score-dup->full_junction_score < config.min_score_diff) {
				filters.push_back("WEAK_SPLIT_ALIGNMENT");
			}
            if ((dup->lc_consensus == NULL || dup->lc_consensus->max_mapq < config.high_confidence_mapq) &&
				(dup->rc_consensus == NULL || dup->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
			}
            if (dup->len() >= config.min_size_for_depth_filtering && dup->disc_pairs < 3) {
                filters.push_back("NOT_ENOUGH_OW_PAIRS");
            }
            if (dup->sr_remap_info && (dup->sr_remap_info->weak_supporting_reads < 3 || dup->sr_remap_info->strong_supporting_reads < 1)) {
				filters.push_back("WEAK_SUPPORT");
			}

            if (filters.empty()) {
                filters.push_back("PASS");
            }

            dup2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, dup, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
        }
    }

    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	auto& bcf_entries_contig = bcf_entries[contig_name];
    	std::sort(bcf_entries_contig.begin(), bcf_entries_contig.end(), [](const bcf1_t* b1, const bcf1_t* b2) {return b1->pos < b2->pos;});
		for (bcf1_t* bcf_entry : bcf_entries_contig) {
 			if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
			}
		}
    }

    chr_seqs.clear();

    flog.close();

    bcf_hdr_destroy(out_vcf_header);
    bcf_close(out_vcf_file);

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
