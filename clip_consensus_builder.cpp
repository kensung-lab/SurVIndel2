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
#include "extend_1sr_consensus.h"

config_t config;
stats_t stats;
std::string workdir, complete_bam_fname, reference_fname;
std::mutex out_mtx, log_mtx;

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

std::unordered_map<std::string, std::vector<deletion_t*> > deletions_by_chr;
std::unordered_map<std::string, std::vector<duplication_t*> > duplications_by_chr;
std::unordered_map<std::string, std::vector<consensus_t*> > unpaired_consensuses_by_chr;
std::vector<double> max_allowed_frac_normalized;
std::mutex mtx, indel_out_mtx, up_consensus_mtx;


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
	indel_t* indel = remap_consensus(joined_consensus, chr_seq, chr_len, ref_lh_start, ref_lh_len, ref_rh_start, ref_rh_len, aligner,
			lc_consensus, rc_consensus, source);
	ref_lh_end = ref_lh_start + ref_lh_len;
	ref_rh_end = ref_rh_start + ref_rh_len;
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


void build_hsr_consensuses(int id, int contig_id, std::string contig_name, std::vector<consensus_t*>& rc_hsr_consensuses,
	std::vector<consensus_t*>& lc_hsr_consensuses) {
    std::string bam_fname = workdir + "/workspace/" + std::to_string(contig_id) + "-HSR.bam";
    open_samFile_t* bam_file = open_samFile(bam_fname, true);

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::deque<bam1_t*> lc_cluster, rc_cluster;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {

        std::pair<int, int> left_and_right_diffs = compute_left_and_right_differences(read, true);

        bool lc_clipped = left_and_right_diffs.first > left_and_right_diffs.second;
        std::deque<bam1_t*>& cluster = lc_clipped ? lc_cluster : rc_cluster;
        std::vector<consensus_t*>& consensuses_v = lc_clipped ? lc_hsr_consensuses : rc_hsr_consensuses;

        if (cluster.size() >= 3 && bam_endpos(cluster.front())-read->core.pos < read->core.l_qseq/2) {
			std::vector<bam_redux_t*> cluster_v;
			for (bam1_t* r : cluster) cluster_v.push_back(new bam_redux_t(r));
			std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, cluster_v, lc_clipped);
			for (auto consensus : consensuses) {
				consensus->is_hsr = true;
				consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
				consensuses_v.push_back(consensus);
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
		std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, cluster_v, false);
		for (auto consensus : consensuses) {
			consensus->is_hsr = true;
			consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
			rc_hsr_consensuses.push_back(consensus);
		}
		for (bam_redux_t* r : cluster_v) delete r;
    }
    for (bam1_t* r : rc_cluster) bam_destroy1(r);
    if (lc_cluster.size() >= 3) {
		std::vector<bam_redux_t*> cluster_v;
		for (bam1_t* r : lc_cluster) cluster_v.push_back(new bam_redux_t(r));
		std::vector<consensus_t*> consensuses = build_full_consensus(contig_id, cluster_v, true);
		for (auto consensus : consensuses) {
			consensus->is_hsr = true;
			consensus->clip_len = consensus_t::UNKNOWN_CLIP_LEN;
			lc_hsr_consensuses.push_back(consensus);
		}
		for (bam_redux_t* r : cluster_v) delete r;
    }
    for (bam1_t* r : lc_cluster) bam_destroy1(r);

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

		for (auto& iv : compatible_lc_idxs) {
			consensus_t* lc_anchor = lc_consensuses[iv.value];
			if (lc_anchor->max_mapq < config.high_confidence_mapq && rc_anchor->max_mapq < config.high_confidence_mapq) continue;

			int min_overlap = (lc_anchor->is_hsr && rc_anchor->is_hsr) ? 50 : std::min(rc_anchor->clip_len, lc_anchor->clip_len)+config.min_clip_len;
			double max_mm_rate = (lc_anchor->is_hsr && rc_anchor->is_hsr) ? 0 : config.max_seq_error;

			suffix_prefix_aln_t spa = aln_suffix_prefix(rc_anchor->consensus, lc_anchor->consensus, 1, -4, min_overlap);
			if (spa.overlap > 0 && spa.mismatch_rate() <= max_mm_rate) {
				rc_lc_scored_pairs.push_back(pair_w_score_t(i, iv.value, spa.score, spa.mismatch_rate(), false));
			} else { // trim low quality (i.e., supported by less than 2 reads) bases
				// TODO: investigate if we can use base qualities for this
				std::string rc_consensus_trim = rc_anchor->consensus.substr(0, rc_anchor->consensus.length()-rc_anchor->lowq_clip_portion);
				std::string lc_consensus_trim = lc_anchor->consensus.substr(lc_anchor->lowq_clip_portion);

				suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus_trim, lc_consensus_trim, 1, -4, min_overlap);
				if (spa.overlap > 0 && spa.mismatch_rate() <= max_mm_rate) {
					rc_lc_scored_pairs.push_back(pair_w_score_t(i, iv.value, spa.score, spa.mismatch_rate(), true));
				}
			}
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
		} else if (indel->indel_type() == "DUP") {
			contig_duplications.push_back((duplication_t*) indel);
		}
	}

	remove_marked_consensuses(rc_consensuses, used_consensus_rc);
	remove_marked_consensuses(lc_consensuses, used_consensus_lc);
}

indel_t* find_indel_from_lc_consensus(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree, std::string contig_name,
		StripedSmithWaterman::Aligner& aligner, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {
	if (!consensus->left_clipped) return NULL;

	hts_pos_t contig_len = chr_seqs.get_len(contig_name);

	hts_pos_t left_ext_target_start = consensus->left_ext_target_start(config);
	hts_pos_t left_ext_target_end = consensus->left_ext_target_end(config);

	hts_pos_t right_ext_target_start = consensus->right_ext_target_start(config);
	hts_pos_t right_ext_target_end = consensus->right_ext_target_end(config);

	int old_consensus_len = consensus->consensus.length();
	extend_consensus_to_left(consensus, candidate_reads_itree, left_ext_target_start, left_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int l_ext_len = consensus->consensus.length() - old_consensus_len;

	old_consensus_len = consensus->consensus.length();
	extend_consensus_to_right(consensus, candidate_reads_itree, right_ext_target_start, right_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int r_ext_len = consensus->consensus.length() - old_consensus_len;

	hts_pos_t remap_target_start = consensus->remap_boundary;
	hts_pos_t remap_target_end = consensus->remap_boundary + config.max_is;
	if (consensus->remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
		remap_target_start = consensus->breakpoint - config.max_is - 2*consensus->consensus.length();
		remap_target_end = consensus->breakpoint + config.max_is + 2*consensus->consensus.length();
	}
	remap_target_start = std::max(remap_target_start-l_ext_len, hts_pos_t(0));
	remap_target_end = std::min(remap_target_end+r_ext_len, contig_len);
	hts_pos_t remap_target_len = remap_target_end - remap_target_start;

	// do not attempt if reference region has Ns - this is because of in our aligner, Ns will always match
	if (remap_target_start >= remap_target_end ||
		has_Ns(chr_seqs.get_seq(contig_name), remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

	if (consensus->clip_len < config.min_clip_len) return NULL;

	hts_pos_t ref_rh_end = consensus->end + EXTRA_SEQ;
	if (ref_rh_end >= contig_len) ref_rh_end = contig_len-1;
	hts_pos_t ref_rh_start = consensus->start - config.max_is;
	if (ref_rh_start < 0) ref_rh_start = 0;
	hts_pos_t ref_rh_len = ref_rh_end - ref_rh_start;

	indel_t* indel = remap_consensus(consensus->consensus, chr_seqs.get_seq(contig_name), contig_len,
			remap_target_start, remap_target_len, ref_rh_start, ref_rh_len, aligner, consensus, NULL, "1SR_LC", false);
	if (indel == NULL) return NULL;

	if (indel->indel_type() == "DUP") {
		duplication_t* dup = (duplication_t*) indel;
		dup->original_start = dup->start, dup->original_end = dup->end;
	}
	indel->extra_info += consensus->name();
	return indel;
}

indel_t* find_indel_from_rc_consensus(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree, std::string contig_name,
		StripedSmithWaterman::Aligner& aligner, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {
	if (consensus->left_clipped) return NULL;

	hts_pos_t contig_len = chr_seqs.get_len(contig_name);

	hts_pos_t left_ext_target_start = consensus->left_ext_target_start(config);
	hts_pos_t left_ext_target_end = consensus->left_ext_target_end(config);

	hts_pos_t right_ext_target_start = consensus->right_ext_target_start(config);
	hts_pos_t right_ext_target_end = consensus->right_ext_target_end(config);

	int old_consensus_len = consensus->consensus.length();
	extend_consensus_to_right(consensus, candidate_reads_itree, right_ext_target_start, right_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int r_ext_len = consensus->consensus.length() - old_consensus_len;

	old_consensus_len = consensus->consensus.length();
	extend_consensus_to_left(consensus, candidate_reads_itree, left_ext_target_start, left_ext_target_end, contig_name, contig_len, config, mateseqs_w_mapq);
	int l_ext_len = consensus->consensus.length() - old_consensus_len;

	hts_pos_t remap_target_end = consensus->remap_boundary;
	hts_pos_t remap_target_start = consensus->remap_boundary - config.max_is;
	if (consensus->remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
		remap_target_start = consensus->breakpoint - config.max_is - 2*consensus->consensus.length();
		remap_target_end = consensus->breakpoint + config.max_is + 2*consensus->consensus.length();
	}
	remap_target_start = std::max(remap_target_start-l_ext_len, hts_pos_t(0));
	remap_target_end = std::min(remap_target_end+r_ext_len, contig_len);
	hts_pos_t remap_target_len = remap_target_end - remap_target_start;

	if (remap_target_start >= remap_target_end ||
		has_Ns(chr_seqs.get_seq(contig_name), remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

	std::string clip = consensus->consensus.substr(consensus->anchor_len());
	if (clip.size() < config.min_clip_len) return NULL;

	hts_pos_t ref_lh_start = consensus->start - EXTRA_SEQ;
	if (ref_lh_start < 0) ref_lh_start = 0;
	hts_pos_t ref_lh_end = consensus->start + consensus->consensus.length();
	if (ref_lh_end >= contig_len) ref_lh_end = contig_len-1;
	hts_pos_t ref_lh_len = ref_lh_end - ref_lh_start;

	indel_t* indel = remap_consensus(consensus->consensus, chr_seqs.get_seq(contig_name), contig_len,
			ref_lh_start, ref_lh_len, remap_target_start, remap_target_len, aligner, NULL,
			consensus, "1SR_RC", true);
	if (indel == NULL) return NULL;

	if (indel->indel_type() == "DUP") {
		duplication_t* dup = (duplication_t*) indel;
		dup->original_start = dup->start, dup->original_end = dup->end;
	}
	indel->extra_info += consensus->name();
	return indel;
}

void find_indels_from_unpaired_consensuses(int id, std::string contig_name, std::vector<consensus_t*>* consensuses,
		int start, int end, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq) {

	out_mtx.lock();
	std::cout << "Dealing with unpaired consensuses in " << contig_name << ", " << start << " to " << end << std::endl;
	out_mtx.unlock();

	std::vector<deletion_t*> local_dels;
	std::vector<duplication_t*> local_dups;

	hts_pos_t contig_len = chr_seqs.get_len(contig_name);

	StripedSmithWaterman::Aligner aligner(1,4,6,1,true);
	open_samFile_t* bam_file = open_samFile(complete_bam_fname);

	std::vector<hts_pair_pos_t> target_ivals;
	for (int i = start; i < consensuses->size() && i < end; i++) {
		consensus_t* consensus = consensuses->at(i);
		hts_pos_t left_ext_target_start = consensus->left_ext_target_start(config);
		hts_pos_t left_ext_target_end = consensus->left_ext_target_end(config);
		hts_pos_t right_ext_target_start = consensus->right_ext_target_start(config);
		hts_pos_t right_ext_target_end = consensus->right_ext_target_end(config);
		hts_pair_pos_t target_ival;
		target_ival.beg = left_ext_target_start, target_ival.end = left_ext_target_end;
		target_ivals.push_back(target_ival);
		target_ival.beg = right_ext_target_start, target_ival.end = right_ext_target_end;
		target_ivals.push_back(target_ival);
	}
	std::vector<ext_read_t*> candidate_reads_for_extension = get_extension_reads(contig_name, target_ivals, contig_len, *mateseqs_w_mapq, config, bam_file);

	std::vector<Interval<ext_read_t*>> it_ivals;
	for (ext_read_t* ext_read : candidate_reads_for_extension) {
		Interval<ext_read_t*> it_ival(ext_read->start, ext_read->end, ext_read);
		it_ivals.push_back(it_ival);
	}
	IntervalTree<ext_read_t*> candidate_reads_for_extension_itree(it_ivals);

	for (int i = start; i < consensuses->size() && i < end; i++) {
		consensus_t* consensus = consensuses->at(i);

		indel_t* smallest_indel = NULL;
		if (!consensus->left_clipped && !consensus->is_hsr) {
			smallest_indel = find_indel_from_rc_consensus(consensus, candidate_reads_for_extension_itree, contig_name, aligner, *mateseqs_w_mapq);
		} else if (consensus->left_clipped && !consensus->is_hsr) {
			smallest_indel = find_indel_from_lc_consensus(consensus, candidate_reads_for_extension_itree, contig_name, aligner, *mateseqs_w_mapq);
		} else if (!consensus->left_clipped && consensus->is_hsr) {
			smallest_indel = remap_rc_cluster(consensus, candidate_reads_for_extension_itree, contig_name, chr_seqs.get_seq(contig_name), contig_len,
					aligner, *mateseqs_w_mapq);
		} else if (consensus->left_clipped && consensus->is_hsr) {
			smallest_indel = remap_lc_cluster(consensus, candidate_reads_for_extension_itree, contig_name, chr_seqs.get_seq(contig_name), contig_len,
					aligner, *mateseqs_w_mapq);
		}

		if (smallest_indel == NULL) {
			delete consensus;
			continue;
		}
		if (smallest_indel->indel_type() == "DEL") {
			local_dels.push_back((deletion_t*) smallest_indel);
		} else {
			local_dups.push_back((duplication_t*) smallest_indel);
		}
	}
	close_samFile(bam_file);
	for (ext_read_t* read : candidate_reads_for_extension) delete read;

	indel_out_mtx.lock();
	deletions_by_chr[contig_name].insert(deletions_by_chr[contig_name].end(), local_dels.begin(), local_dels.end());
	duplications_by_chr[contig_name].insert(duplications_by_chr[contig_name].end(), local_dups.begin(), local_dups.end());
	indel_out_mtx.unlock();
}

void build_clip_consensuses(int id, int contig_id, std::string contig_name, std::vector<deletion_t*>& deletions,
                            std::vector<duplication_t*>& duplications, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	std::vector<deletion_t*> contig_deletions;
	std::vector<duplication_t*> contig_duplications;

	StripedSmithWaterman::Aligner aligner(1,4,6,1,true);
	StripedSmithWaterman::Filter filter;

    std::string clip_fname = workdir + "/workspace/" + std::to_string(contig_id) + "-CLIP.bam";
    open_samFile_t* clip_file = open_samFile(clip_fname, true);

    hts_itr_t* iter = sam_itr_querys(clip_file->idx, clip_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::vector<bam_redux_t*> lc_reads, rc_reads;
    while (sam_itr_next(clip_file->file, iter, read) >= 0) {
        if (is_left_clipped(read, config.min_clip_len)) {
            lc_reads.push_back(new bam_redux_t(read));
        }
        if (is_right_clipped(read, config.min_clip_len)) {
            rc_reads.push_back(new bam_redux_t(read));
        }
    }

    bam_destroy1(read);

    hts_itr_destroy(iter);
    close_samFile(clip_file);

    std::vector<consensus_t*> rc_consensuses, lc_consensuses;
    std::vector<bam_redux_t*> curr_candidate_cluster;

    auto lc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->start-r2->start) <= config.max_clipped_pos_dist;};
    if (!lc_reads.empty()) {
        for (bam_redux_t* lc_read : lc_reads) {
            if (!curr_candidate_cluster.empty() &&
                !lc_same_cluster(curr_candidate_cluster[0], lc_read)) { // candidate cluster complete
                std::vector<consensus_t*> full_consensuses = build_full_consensus(contig_id, curr_candidate_cluster, true);
				lc_consensuses.insert(lc_consensuses.end(), full_consensuses.begin(), full_consensuses.end());
                curr_candidate_cluster.clear();
            }
            curr_candidate_cluster.push_back(lc_read);
        }
        // process last cluster
        std::vector<consensus_t*> full_consensuses = build_full_consensus(contig_id, curr_candidate_cluster, true);
        curr_candidate_cluster.clear();
        lc_consensuses.insert(lc_consensuses.end(), full_consensuses.begin(), full_consensuses.end());
    }
    auto rc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->end-r2->end) <= config.max_clipped_pos_dist;};
    sort(rc_reads.begin(), rc_reads.end(), [](bam_redux_t* r1, bam_redux_t* r2) { return r1->end < r2->end; });
    if (!rc_reads.empty()) {
        for (bam_redux_t* rc_read : rc_reads) {
            if (!curr_candidate_cluster.empty() &&
                !rc_same_cluster(curr_candidate_cluster[0], rc_read)) { // candidate cluster complete
            	std::vector<consensus_t*> full_consensuses = build_full_consensus(contig_id, curr_candidate_cluster, false);
				rc_consensuses.insert(rc_consensuses.end(), full_consensuses.begin(), full_consensuses.end());
                curr_candidate_cluster.clear();
            }
            curr_candidate_cluster.push_back(rc_read);
        }
        // process last cluster
        std::vector<consensus_t*> full_consensuses = build_full_consensus(contig_id, curr_candidate_cluster, false);
        curr_candidate_cluster.clear();
		rc_consensuses.insert(rc_consensuses.end(), full_consensuses.begin(), full_consensuses.end());
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

	// find RC-LC pairs
	find_indels_from_rc_lc_pairs(contig_name, rc_consensuses, lc_consensuses, contig_deletions, contig_duplications, aligner);
	find_indels_from_rc_lc_pairs(contig_name, rc_hsr_consensuses, lc_consensuses, contig_deletions, contig_duplications, aligner);
	find_indels_from_rc_lc_pairs(contig_name, rc_consensuses, lc_hsr_consensuses, contig_deletions, contig_duplications, aligner);
	find_indels_from_rc_lc_pairs(contig_name, rc_hsr_consensuses, lc_hsr_consensuses, contig_deletions, contig_duplications, aligner);

	deletions.insert(deletions.end(), contig_deletions.begin(), contig_deletions.end());
	duplications.insert(duplications.end(), contig_duplications.begin(), contig_duplications.end());

    /* == Deal with unpaired consensuses == */
	up_consensus_mtx.lock();
	std::vector<consensus_t*>& unpaired_consensuses = unpaired_consensuses_by_chr[contig_name];
	unpaired_consensuses.insert(unpaired_consensuses.end(), rc_consensuses.begin(), rc_consensuses.end());
	unpaired_consensuses.insert(unpaired_consensuses.end(), lc_consensuses.begin(), lc_consensuses.end());
	unpaired_consensuses.insert(unpaired_consensuses.end(), rc_hsr_consensuses.begin(), rc_hsr_consensuses.end());
	unpaired_consensuses.insert(unpaired_consensuses.end(), lc_hsr_consensuses.begin(), lc_hsr_consensuses.end());
	up_consensus_mtx.unlock();
}


void build_consensuses(int id, int contig_id, std::string contig_name, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq) {
    mtx.lock();
    std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
    std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
    mtx.unlock();

    mtx.lock();
	std::cout << "Computing indels for " << contig_name << std::endl;
	mtx.unlock();

	build_clip_consensuses(id, contig_id, contig_name, deletions, duplications, *mateseqs_w_mapq);
}

void size_and_depth_filtering(int id, std::string contig_name) {
    open_samFile_t* bam_file = get_bam_reader(complete_bam_fname);

    out_mtx.lock();
    std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
    out_mtx.unlock();
    std::vector<double> temp1;
    std::vector<uint32_t> temp2;
    calculate_confidence_interval_size(contig_name, temp1, temp2, deletions, bam_file, config, stats);
    depth_filter_del(contig_name, deletions, bam_file, config.min_size_for_depth_filtering, config);

    out_mtx.lock();
    std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
    out_mtx.unlock();
//    std::vector<duplication_t*> duplications_w_cleanup, duplications_wo_cleanup;
//    for (duplication_t* dup : duplications) {
//        if (dup->len() >= config.min_size_for_depth_filtering) duplications_w_cleanup.push_back(dup);
//        else duplications_wo_cleanup.push_back(dup);
//    }
//    depth_filter_dup_w_cleanup(contig_name, duplications_w_cleanup, bam_file, stats, config, max_allowed_frac_normalized, workdir);
//    depth_filter_dup(contig_name, duplications_wo_cleanup, bam_file, config.min_size_for_depth_filtering, config);
	depth_filter_dup(contig_name, duplications, bam_file, config.min_size_for_depth_filtering, config);
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

    if (config.log) flog.open(workdir + "/log");
    else flog.open("/dev/null");

    std::vector<std::unordered_map<std::string, std::pair<std::string, int> > > mateseqs_w_mapq(contig_map.size());

    ctpl::thread_pool thread_pool1(config.threads);
    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string fname = workdir + "/workspace/" + std::to_string(contig_id) + ".dc_mateseqs";
		std::ifstream fin(fname);
		std::string qname, read_seq; int mapq;
		while (fin >> qname >> read_seq >> mapq) {
			mateseqs_w_mapq[contig_id][qname] = {read_seq, mapq};
		}

		std::string contig_name = contig_map.get_name(contig_id);
//		if (contig_name != "chr22") continue;
        std::future<void> future = thread_pool1.push(build_consensuses, contig_id, contig_name, &mateseqs_w_mapq[contig_id]);
        futures.push_back(std::move(future));
    }
    thread_pool1.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

    int block_size = 1000;
    ctpl::thread_pool thread_pool2(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
    	std::string contig_name = contig_map.get_name(contig_id);
    	std::vector<consensus_t*>& consensuses = unpaired_consensuses_by_chr[contig_name];
    	if (consensuses.empty()) continue;
    	for (int i = 0; i <= consensuses.size()/block_size; i++) {
    		int start = i * block_size;
    		int end = std::min(start+block_size, (int) consensuses.size());
			std::future<void> future = thread_pool2.push(find_indels_from_unpaired_consensuses,
					contig_name, &consensuses, i*block_size, end, &mateseqs_w_mapq[contig_id]);
			futures.push_back(std::move(future));
    	}
    }
    thread_pool2.stop(true);
	for (size_t i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		auto& deletions = deletions_by_chr[contig_name];
		for (int i = 0; i < deletions.size(); i++) {
			if (deletions[i]->len() < config.min_sv_size) {
				delete deletions[i];
				deletions[i] = NULL;
			}
		}
		deletions_by_chr[contig_name].erase(std::remove(deletions_by_chr[contig_name].begin(), deletions_by_chr[contig_name].end(), (deletion_t*) NULL), deletions_by_chr[contig_name].end());
		auto& duplications = duplications_by_chr[contig_name];
		for (int i = 0; i < duplications.size(); i++) {
			if (duplications[i]->len() < config.min_sv_size) {
				delete duplications[i];
				duplications[i] = NULL;
			}
		}
		duplications_by_chr[contig_name].erase(std::remove(duplications_by_chr[contig_name].begin(), duplications_by_chr[contig_name].end(), (duplication_t*) NULL), duplications_by_chr[contig_name].end());
	}

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

    ctpl::thread_pool thread_pool3(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool3.push(size_and_depth_filtering, contig_name);
        futures.push_back(std::move(future));
    }
    thread_pool3.stop(true);
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
            if (del->full_junction_score && del->lh_best1_junction_score+del->rh_best1_junction_score-del->full_junction_score < config.min_score_diff) {
            	filters.push_back("WEAK_SPLIT_ALIGNMENT");
            }
            if ((del->lc_consensus == NULL || del->lc_consensus->max_mapq < config.high_confidence_mapq) &&
            	(del->rc_consensus == NULL || del->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
            }
            if (del->source == "1SR_RC" || del->source == "1HSR_RC") {
            	if (del->rc_consensus->right_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            } else if (del->source == "1SR_LC" || del->source == "1HSR_LC") {
            	if (del->lc_consensus->left_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            }

            if (filters.empty()) {
                filters.push_back("PASS");
            }

            del2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, del, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
            delete del->rc_consensus;
            delete del->lc_consensus;
            delete del;
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
            if (dup->full_junction_score && dup->lh_best1_junction_score+dup->rh_best1_junction_score-dup->full_junction_score < config.min_score_diff) {
				filters.push_back("WEAK_SPLIT_ALIGNMENT");
			}
            if ((dup->lc_consensus == NULL || dup->lc_consensus->max_mapq < config.high_confidence_mapq) &&
				(dup->rc_consensus == NULL || dup->rc_consensus->max_mapq < config.high_confidence_mapq)) {
				filters.push_back("LOW_MAPQ_CONSENSUSES");
			}
            if (dup->len() >= config.min_size_for_depth_filtering && dup->disc_pairs < 3) {
                filters.push_back("NOT_ENOUGH_OW_PAIRS");
            }
            if (dup->source == "1SR_RC" || dup->source == "1HSR_RC") {
				if (dup->rc_consensus->right_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
            } else if (dup->source == "1SR_LC" || dup->source == "1HSR_LC") {
            	if (dup->lc_consensus->left_ext_reads < 3) filters.push_back("FAILED_TO_EXTEND");
			}

            if (filters.empty()) {
                filters.push_back("PASS");
            }

            dup2bcf(out_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, dup, filters);
            bcf_entries[contig_name].push_back(bcf_dup(bcf_entry));
            delete dup->rc_consensus;
            delete dup->lc_consensus;
            delete dup;
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
