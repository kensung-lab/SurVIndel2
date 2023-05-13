#ifndef REMAPPING_H_
#define REMAPPING_H_

#include "utils.h"
#include "libs/ssw_cpp.h"
#include "libs/ssw.h"
#include <mutex>

extern config_t config;

std::vector<StripedSmithWaterman::Alignment> get_best_alns(char* ref, int ref_len, char* query, StripedSmithWaterman::Aligner& aligner) {
	std::vector<StripedSmithWaterman::Alignment> best_alns;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	aligner.Align(query, ref, ref_len, filter, &aln, 15);
	best_alns.push_back(aln);
	int best_score = aln.sw_score;
	int new_start = aln.ref_begin+5;

	if (new_start >= ref_len) return best_alns;

	aligner.Align(query, ref+new_start, ref_len-new_start, filter, &aln, 15);
	while (aln.sw_score == best_score) {
		aln.ref_begin += new_start;
		aln.ref_end += new_start;
		best_alns.push_back(aln);
		new_start += aln.ref_begin+5;
		if (new_start >= ref_len) break;
		aligner.Align(query, ref+new_start, ref_len-new_start, filter, &aln, 15);
	}
	return best_alns;
}

indel_t* remap_consensus(std::string& consensus_seq, char* reference, int reference_len, hts_pos_t& ref_lh_start, hts_pos_t& ref_lh_len,
		hts_pos_t& ref_rh_start, hts_pos_t& ref_rh_len, StripedSmithWaterman::Aligner& aligner, consensus_t* lc_consensus,
		consensus_t* rc_consensus, std::string source, bool choose_rightmost_indel = false) {

	if (has_Ns(reference, ref_lh_start, ref_lh_len) || has_Ns(reference, ref_rh_start, ref_rh_len)) return NULL;

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment lh_full_aln, rh_full_aln;
	aligner.Align(consensus_seq.c_str(), reference+ref_lh_start, ref_lh_len, filter, &lh_full_aln, 0);
	aligner.Align(consensus_seq.c_str(), reference+ref_rh_start, ref_rh_len, filter, &rh_full_aln, 0);
	StripedSmithWaterman::Alignment& best_aln = lh_full_aln.sw_score >= rh_full_aln.sw_score ? lh_full_aln : rh_full_aln;

	char ref_lh_cstr[100000];
	for (int i = 0; i < ref_lh_len; i++) {
		ref_lh_cstr[i] = toupper(reference[ref_lh_start+i]);
	} ref_lh_cstr[ref_lh_len] = '\0';
	int* prefix_scores = smith_waterman_gotoh(ref_lh_cstr, ref_lh_len, consensus_seq.c_str(), consensus_seq.length(), 1, -4, -6, -1);

	char ref_rh_cstr[100000];
	for (int i = 0; i < ref_rh_len; i++) {
		ref_rh_cstr[i] = toupper(reference[ref_rh_start+ref_rh_len-1-i]);
	} ref_rh_cstr[ref_rh_len] = '\0';
	std::string consensus_seq_rev = std::string(consensus_seq.rbegin(), consensus_seq.rend());
	int* suffix_scores = smith_waterman_gotoh(ref_rh_cstr, ref_rh_len, consensus_seq_rev.c_str(), consensus_seq_rev.length(), 1, -4, -6, -1);

	int max_score = 0, split_i = 0;
	for (int i = config.min_clip_len; i < consensus_seq.length()-config.min_clip_len; i++) {
		int prefix_score = prefix_scores[i-1], suffix_score = suffix_scores[consensus_seq.length()-i-1];
		if (prefix_score + suffix_score > max_score || (prefix_score + suffix_score == max_score && choose_rightmost_indel)) {
			max_score = prefix_score + suffix_score;
			split_i = i;
		}
	}

	delete[] prefix_scores;
	delete[] suffix_scores;

	char* consensus_cstr = new char[consensus_seq.length()+1];
	strcpy(consensus_cstr, consensus_seq.c_str());
	char terminator = '\0';
	std::swap(terminator, consensus_cstr[split_i]);

	StripedSmithWaterman::Alignment lh_aln, rh_aln;
	aligner.Align(consensus_cstr, reference+ref_lh_start, ref_lh_len, filter, &lh_aln, 0);
	if (lh_aln.ref_begin - get_left_clip_size(lh_aln) < 0) {
		int extend_by = std::min((int) ref_lh_start, get_left_clip_size(lh_aln) - lh_aln.ref_begin);
		ref_lh_start -= extend_by;
		ref_lh_len += extend_by;
	}
	std::vector<StripedSmithWaterman::Alignment> lh_alns = get_best_alns(reference+ref_lh_start, ref_lh_len, consensus_cstr, aligner);
	std::swap(terminator, consensus_cstr[split_i]);

	// determine if we need to extend the area
	aligner.Align(consensus_cstr+split_i, reference+ref_rh_start, ref_rh_len, filter, &rh_aln, 0);
	if (rh_aln.ref_end + get_right_clip_size(rh_aln) >= ref_rh_len) {
		ref_rh_len = std::min(rh_aln.ref_end + get_right_clip_size(rh_aln), int(reference_len-ref_rh_start));
	}
	std::vector<StripedSmithWaterman::Alignment> rh_alns = get_best_alns(reference+ref_rh_start, ref_rh_len, consensus_cstr+split_i, aligner);
	delete[] consensus_cstr;

	int min_size = INT32_MAX;
	for (StripedSmithWaterman::Alignment& _lh_aln : lh_alns) {
		for (StripedSmithWaterman::Alignment& _rh_aln : rh_alns) {
			int size = abs((ref_rh_start + _rh_aln.ref_begin - 1) - (ref_lh_start + _lh_aln.ref_end));
			if (min_size > size || (min_size == size && choose_rightmost_indel)) {
				lh_aln = _lh_aln, rh_aln = _rh_aln;
				min_size = size;
			}
		}
	}

	hts_pos_t la_start = ref_lh_start + lh_aln.ref_begin, ra_end = ref_rh_start + rh_aln.ref_end;
	int start = ref_lh_start + lh_aln.ref_end, end = ref_rh_start + rh_aln.ref_begin - 1;

	if (start == end) return NULL;
	indel_t* indel;
	if (start < end) {
		std::string ins_seq = consensus_seq.substr(lh_aln.query_end, split_i-lh_aln.query_end-1);
		ins_seq += consensus_seq.substr(split_i, rh_aln.query_begin);
		indel = new deletion_t(start, end, la_start, ra_end, lc_consensus, rc_consensus, lh_aln.sw_score, rh_aln.sw_score, source, ins_seq);
	} else {
		int overlap_len = start - end;
		int lh_suffix_score = find_aln_suffix_score(lh_aln.cigar, overlap_len, 1, -4, -6, -1);
		int rh_prefix_score = find_aln_prefix_score(rh_aln.cigar, overlap_len, 1, -4, -6, -1);
		if (ra_end-start < config.min_clip_len || end-la_start < config.min_clip_len ||
				(lh_suffix_score == overlap_len && rh_prefix_score == overlap_len)) { // perfect dup
			std::string ins_seq = consensus_seq.substr(split_i, get_left_clip_size(rh_aln));
			if (ins_seq.length() > start-end) return NULL;
			indel = new duplication_t(end, start, la_start, ra_end, lc_consensus, rc_consensus, source, ins_seq);
		} else { // insertion
			if (lh_suffix_score <= rh_prefix_score) {
				std::string ins_seq = consensus_seq.substr(split_i-overlap_len, overlap_len);
				ins_seq += consensus_seq.substr(split_i, get_left_clip_size(rh_aln));
				indel = new duplication_t(end, end, la_start, ra_end, lc_consensus, rc_consensus, source, ins_seq);
			} else {
				int ins_seq_len = overlap_len + get_left_clip_size(rh_aln);
				std::string ins_seq = consensus_seq.substr(split_i, ins_seq_len);
				indel = new duplication_t(start, start, la_start, ra_end, lc_consensus, rc_consensus, source, ins_seq);
			}
		}
	}
	indel->full_junction_score = std::max(lh_full_aln.sw_score, rh_full_aln.sw_score);
	indel->lh_best1_junction_score = lh_aln.sw_score, indel->rh_best1_junction_score = rh_aln.sw_score;
	indel->lh_best2_junction_score = lh_aln.sw_score_next_best, indel->rh_best2_junction_score = rh_aln.sw_score_next_best;
	if (rc_consensus && rc_consensus->clip_len == consensus_t::UNKNOWN_CLIP_LEN) {
		rc_consensus->clip_len = consensus_seq.length() - split_i;
	}
	if (lc_consensus && lc_consensus->clip_len == consensus_t::UNKNOWN_CLIP_LEN) {
		lc_consensus->clip_len = split_i;
	}

	indel->extra_info += lh_aln.cigar_string + "," + rh_aln.cigar_string + "," + best_aln.cigar_string + ",";
	return indel;
}


#endif /* REMAPPING_H_ */
