#ifndef HSR_H_
#define HSR_H_

#include <htslib/sam.h>

#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "utils.h"
#include "sam_utils.h"
#include "remapping.h"
#include "extend_1sr_consensus.h"

extern config_t config;

const int EXTRA_SEQ = 10;


indel_t* remap_rc_cluster(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree, std::string contig_name,
		char* contig_seq, hts_pos_t contig_len, StripedSmithWaterman::Aligner& aligner,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

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

	hts_pos_t ref_start = std::max(hts_pos_t(0), consensus->start - EXTRA_SEQ);
	hts_pos_t ref_end = std::min(consensus->end + config.max_is, contig_len);
	hts_pos_t ref_len = ref_end - ref_start;

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
		has_Ns(contig_seq, remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

    indel_t* indel = remap_consensus(consensus->consensus, contig_seq, contig_len, ref_start, ref_len, remap_target_start,
    		remap_target_len, aligner, NULL, consensus, "1HSR_RC");
    if (indel != NULL) indel->extra_info += consensus->name() + ",";
    return indel;
}

indel_t* remap_lc_cluster(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree, std::string contig_name,
		char* contig_seq, hts_pos_t contig_len, StripedSmithWaterman::Aligner& aligner,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

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

	hts_pos_t ref_end = std::min(consensus->end + EXTRA_SEQ, contig_len);
	hts_pos_t ref_start = std::max(hts_pos_t(0), consensus->start - config.max_is);
	hts_pos_t ref_len = ref_end - ref_start;

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
		has_Ns(contig_seq, remap_target_start, remap_target_end-remap_target_start)) {
		return NULL;
	}

	indel_t* indel = remap_consensus(consensus->consensus, contig_seq, contig_len, remap_target_start, remap_target_len,
			ref_start, ref_len, aligner, consensus, NULL, "1HSR_LC");
	if (indel != NULL) indel->extra_info += consensus->name() + ",";
    return indel;
}

// remove HSR clusters that overlap with a clipped position, i.e., the clipped position of a clipped cluster is contained in the HSR cluster
// direction of clip must be the same
// clipped_consensuses must be sorted by breakpoint
void remove_hsr_overlapping_clipped(std::vector<consensus_t*>& hsr_consensuses, std::vector<consensus_t*>& clipped_consensuses) {
	std::sort(hsr_consensuses.begin(), hsr_consensuses.end(),
			[](const consensus_t* c1, const consensus_t* c2) { return c1->start < c2->start; });
	std::vector<consensus_t*> kept_consensus;
	int i = 0;
	for (consensus_t* c : hsr_consensuses) {
		while (i < clipped_consensuses.size() && clipped_consensuses[i]->breakpoint < c->start) i++;

		// clipped_consensuses[i] must be same clip direction as cluster and breakpoint must be the smallest s.t. >= cluster.start
		if (i >= clipped_consensuses.size() || clipped_consensuses[i]->breakpoint > c->end) {
			kept_consensus.push_back(c);
		}
	}
	kept_consensus.swap(hsr_consensuses);
}

int compute_left_half_Ms(bam1_t* r) {
	int border = r->core.l_qseq/2;
	int left_Ms = 0;

	// compute how many Ms (either matches or mismatches are in the first half of the read)
	// in other words, this is readlen/2 - number of insertions in the first half
	int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
	for (int i = 0; i < r->core.n_cigar; i++) {
		char op_char = bam_cigar_opchr(cigar[i]);
		int op_len = bam_cigar_oplen(cigar[i]);
		if (op_char == 'M') {
			left_Ms += std::min(border, qpos+op_len) - qpos;
			qpos += op_len;
		} else if (op_char == 'I') {
			qpos += op_len;
		}

		if (qpos >= border) break;
	}
	return left_Ms;
}

std::pair<int, int> compute_left_and_right_differences_indel_as_1_diff(bam1_t* r) {
	int border = r->core.l_qseq/2;
	int left_Ms = compute_left_half_Ms(r);

	// we computed left_Ms because the MD tag does not include insertions
	// i.e. if I have CIGAR: 50M50I50M and MD: 74A25, the mismatch is not in position 75 in the query,
	// but it is in position 125
	// here we compute the number of deleted bases and mismatches in each half of the read based on the MD tag
	std::string md_tag = bam_aux2Z(bam_aux_get(r, "MD"));
	int m = 0;
	int left_diffs = 0, right_diffs = 0;
	bool del_mode = false;
	for (int i = 0; i < md_tag.length(); i++) {
		char c = md_tag[i];
		if (c >= '0' && c <= '9') {
			m = m*10 + c-'0';
			del_mode = false;
		} else if (c >= 'A' && c <= 'T') {
			left_Ms -= m;
			m = 0;
			if (!del_mode) {
				left_Ms--;
				if (left_Ms > 0) left_diffs++;
				else right_diffs++;
			}
		} else if (c == '^') {
			del_mode = true;
			if (left_Ms > 0) left_diffs++;
			else right_diffs++;
		}
	}

	int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
	for (int i = 0; i < r->core.n_cigar; i++) {
		char op_char = bam_cigar_opchr(cigar[i]);
		int op_len = bam_cigar_oplen(cigar[i]);
		if (op_char == 'M') {
			qpos += op_len;
		} else if (op_char == 'I') {
			if (qpos < border && qpos + op_len > border) {  // this is the case there an insertion is partially in the left
															// half and partially in the right half
				left_diffs++;
				right_diffs++;
			} else if (qpos < border) {
				left_diffs++;
			} else if (qpos >= border) {
				right_diffs++;
			}
			qpos += op_len;
		}
	}

	return {left_diffs, right_diffs};
}

// computes the differences for the left and the right half of the read
std::pair<int, int> compute_left_and_right_differences_indel_as_n_diffs(bam1_t* r) {
    int border = r->core.l_qseq/2;
    int left_Ms = compute_left_half_Ms(r);

    // we computed left_Ms because the MD tag does not include insertions
    // i.e. if I have CIGAR: 50M50I50M and MD: 74A25, the mismatch is not in position 75 in the query,
    // but it is in position 125
    // here we compute the number of deleted bases and mismatches in each half of the read based on the MD tag
    std::string md_tag = bam_aux2Z(bam_aux_get(r, "MD"));
    int m = 0;
    int left_diffs = 0, right_diffs = 0;
    bool del_mode = false;
    for (int i = 0; i < md_tag.length(); i++) {
        char c = md_tag[i];
        if (c >= '0' && c <= '9') {
            m = m*10 + c-'0';
            del_mode = false;
        } else if (c >= 'A' && c <= 'T') {
            left_Ms -= m;
            m = 0;
            if (left_Ms > 0) left_diffs++;
            else right_diffs++;
            if (!del_mode) left_Ms--;
        } else if (c == '^') {
            del_mode = true;
        }
    }

    // We also need to count the number of inserted bases, since the MD tag does not report this information
    int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
    for (int i = 0; i < r->core.n_cigar; i++) {
        char op_char = bam_cigar_opchr(cigar[i]);
        int op_len = bam_cigar_oplen(cigar[i]);
        if (op_char == 'M') {
            qpos += op_len;
        } else if (op_char == 'I') {
            if (qpos < border && qpos + op_len > border) {  // this is the case there an insertion is partially in the left
                                                            // half and partially in the right half
                left_diffs += border - qpos;
                right_diffs += qpos + op_len - border;
            } else if (qpos < border) {
                left_diffs += op_len;
            } else if (qpos >= border) {
                right_diffs += op_len;
            }
            qpos += op_len;
        }
    }

    return {left_diffs, right_diffs};
}

std::pair<int, int> compute_left_and_right_differences(bam1_t* r, bool indel_as_single_diff) {
	if (indel_as_single_diff) {
		return compute_left_and_right_differences_indel_as_1_diff(r);
	} else {
		return compute_left_and_right_differences_indel_as_n_diffs(r);
	}
}

void filter_well_aligned_to_ref(char* contig_seq, hts_pos_t contig_len, std::vector<consensus_t*>& consensuses,
                                StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
    std::vector<consensus_t*> retained;

    StripedSmithWaterman::Alignment aln;
    for (consensus_t* c : consensuses) {
        std::string padded_consensus = config.clip_penalty_padding() + c->consensus + config.clip_penalty_padding();
        hts_pos_t ref_start = c->start - 2*config.clip_penalty;
        if (ref_start < 0) ref_start = 0;
        hts_pos_t ref_end = c->end + 2*config.clip_penalty;
        if (ref_end >= contig_len) ref_end = contig_len-1;
        aligner.Align(padded_consensus.c_str(), contig_seq+ref_start+1, ref_end-ref_start, filter, &aln, 0);
        int differences = std::count(aln.cigar_string.begin(), aln.cigar_string.end(), 'D') +
                          std::count(aln.cigar_string.begin(), aln.cigar_string.end(), 'I') + aln.mismatches - 2*config.clip_penalty;
        if ((aln.sw_score < padded_consensus.length()*2 && differences >= config.min_diff_hsr) || is_clipped(aln)) {
            retained.push_back(c);
        } else {
        	delete c;
        }
    }

    consensuses.swap(retained);
}

void select_nonoverlapping_clusters(std::vector<consensus_t*>& consensuses) {
	// for overlapping pairs, keep the ones with higher count
	std::vector<consensus_t*> to_be_deleted;
	std::vector<Interval<consensus_t*> > rc_iv;
	for (consensus_t* c : consensuses) rc_iv.push_back(Interval<consensus_t*>(c->start, c->end, c));
	IntervalTree<consensus_t*> rc_it(rc_iv);
	std::vector<consensus_t*> kept_consensuses;
	for (consensus_t* c : consensuses) {
		std::vector<Interval<consensus_t*> > ov = rc_it.findOverlapping(c->start, c->end);
		bool keep_consensus = true;
		for (Interval<consensus_t*> ov_c : ov) {
			if (ov_c.value->supp_clipped_reads > c->supp_clipped_reads) { // a higher count was found
				keep_consensus = false;
				break;
			}
		}

		if (keep_consensus) {
			kept_consensuses.push_back(c);
		} else {
			to_be_deleted.push_back(c);
		}
	}
	for (consensus_t* c : to_be_deleted) delete c;
	consensuses.swap(kept_consensuses);
}

#endif /* HSR_H_ */
