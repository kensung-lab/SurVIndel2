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
	aligner.Align(query, ref, ref_len, filter, &aln, 0);
	best_alns.push_back(aln);
	int best_score = aln.sw_score;
	int new_start = aln.ref_begin+5;

	aligner.Align(query, ref+new_start, ref_len-new_start, filter, &aln, 0);
	while (aln.sw_score == best_score) {
		aln.ref_begin += new_start;
		aln.ref_end += new_start;
		best_alns.push_back(aln);
		new_start += aln.ref_begin+5;
		if (new_start >= ref_len) break;
		aligner.Align(query, ref+new_start, ref_len-new_start, filter, &aln, 0);
	}
	return best_alns;
}

indel_t* remap_consensus(std::string& consensus_seq, char* reference, hts_pos_t ref_lh_start, hts_pos_t ref_lh_len,
		hts_pos_t ref_rh_start, hts_pos_t ref_rh_len, StripedSmithWaterman::Aligner& aligner, consensus_t* lc_consensus,
		consensus_t* rc_consensus, std::string source, hts_pos_t& la_start, hts_pos_t& ra_end, bool choose_rightmost_indel = false) {

	if (has_Ns(reference, ref_lh_start, ref_lh_len) || has_Ns(reference, ref_rh_start, ref_rh_len)) return NULL;

	std::string padded_consensus = config.clip_penalty_padding() + consensus_seq + config.clip_penalty_padding();

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment lh_full_aln, rh_full_aln;
	aligner.Align(padded_consensus.c_str(), reference+ref_lh_start, ref_lh_len, filter, &lh_full_aln, 0);
	aligner.Align(padded_consensus.c_str(), reference+ref_rh_start, ref_rh_len, filter, &rh_full_aln, 0);
	StripedSmithWaterman::Alignment& best_aln = lh_full_aln.sw_score >= rh_full_aln.sw_score ? lh_full_aln : rh_full_aln;

	char ref_lh_cstr[100000];
	for (int i = 0; i < ref_lh_len; i++) {
		ref_lh_cstr[i] = toupper(reference[ref_lh_start+i]);
	} ref_lh_cstr[ref_lh_len] = '\0';
	int* prefix_scores = smith_waterman_gotoh(ref_lh_cstr, ref_lh_len, padded_consensus.c_str(), padded_consensus.length(), 1, -4, -6, -1);

	char ref_rh_cstr[100000];
	for (int i = 0; i < ref_rh_len; i++) {
		ref_rh_cstr[i] = toupper(reference[ref_rh_start+ref_rh_len-1-i]);
	} ref_rh_cstr[ref_rh_len] = '\0';
	std::string padded_consensus_rev = std::string(padded_consensus.rbegin(), padded_consensus.rend());
	int* suffix_scores = smith_waterman_gotoh(ref_rh_cstr, ref_rh_len, padded_consensus_rev.c_str(), padded_consensus_rev.length(), 1, -4, -6, -1);

	int max_score = 0, split_i = 0;
	for (int i = config.min_clip_len; i < consensus_seq.length()-config.min_clip_len; i++) {
		int prefix_score = prefix_scores[config.clip_penalty+i-1], suffix_score = suffix_scores[consensus_seq.length()-i-1+config.clip_penalty];
		if (prefix_score + suffix_score > max_score || (prefix_score + suffix_score == max_score && choose_rightmost_indel)) {
			max_score = prefix_score + suffix_score;
			split_i = i;
		}
	}

	delete[] prefix_scores;
	delete[] suffix_scores;

	char* padded_consensus_cstr = new char[padded_consensus.length()+1];
	strcpy(padded_consensus_cstr, padded_consensus.c_str());
	char terminator = '\0';
	std::swap(terminator, padded_consensus_cstr[config.clip_penalty+split_i]);

	std::vector<StripedSmithWaterman::Alignment> lh_alns = get_best_alns(reference+ref_lh_start, ref_lh_len, padded_consensus_cstr, aligner);
	std::swap(terminator, padded_consensus_cstr[config.clip_penalty+split_i]);

	std::vector<StripedSmithWaterman::Alignment> rh_alns = get_best_alns(reference+ref_rh_start, ref_rh_len, padded_consensus_cstr+config.clip_penalty+split_i, aligner);
	delete[] padded_consensus_cstr;

	StripedSmithWaterman::Alignment lh_aln, rh_aln;
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

	la_start = ref_lh_start + lh_aln.ref_begin, ra_end = ref_rh_start + rh_aln.ref_end;
	int start = ref_lh_start + lh_aln.ref_end, end = ref_rh_start + rh_aln.ref_begin - 1;

	if (start == end) return NULL;
	indel_t* indel;
	if (start < end) {
		std::string ins_seq = padded_consensus.substr(lh_aln.query_end, config.clip_penalty+split_i-lh_aln.query_end-1);
		ins_seq += padded_consensus.substr(config.clip_penalty+split_i, rh_aln.query_begin);
		indel = new deletion_t(start, end, la_start, ra_end, lc_consensus, rc_consensus, lh_aln.sw_score, rh_aln.sw_score, source, ins_seq);
	} else {
		std::string ins_seq = padded_consensus.substr(config.clip_penalty+split_i, get_left_clip_size(rh_aln));
		if (ins_seq.length() > start-end) return NULL;
		indel = new duplication_t(end, start, la_start, ra_end, lc_consensus, rc_consensus, source, ins_seq);
	}
	indel->full_junction_score = std::max(lh_full_aln.sw_score, rh_full_aln.sw_score);
	indel->split_junction_score = lh_aln.sw_score + rh_aln.sw_score;
	if (rc_consensus && rc_consensus->clip_len == consensus_t::UNKNOWN_CLIP_LEN) {
		rc_consensus->clip_len = consensus_seq.length() - split_i;
	}
	if (lc_consensus && lc_consensus->clip_len == consensus_t::UNKNOWN_CLIP_LEN) {
		lc_consensus->clip_len = split_i;
	}

	indel->extra_info += lh_aln.cigar_string + "," + rh_aln.cigar_string + "," + best_aln.cigar_string + ",";
	return indel;
}

indel_t* remap_consensus(std::string& consensus_seq, char* reference, hts_pos_t ref_lh_start, hts_pos_t ref_lh_len,
		hts_pos_t ref_rh_start, hts_pos_t ref_rh_len, StripedSmithWaterman::Aligner& aligner, consensus_t* lc_consensus,
		consensus_t* rc_consensus, std::string source, bool choose_rightmost_indel = false) {
	hts_pos_t la_start, ra_end;
	return remap_consensus(consensus_seq, reference, ref_lh_start, ref_lh_len, ref_rh_start, ref_rh_len, aligner,
			lc_consensus, rc_consensus, source, la_start, ra_end, choose_rightmost_indel);
}


#endif /* REMAPPING_H_ */
