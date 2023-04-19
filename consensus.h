#ifndef CONSENSUS_H_
#define CONSENSUS_H_

#include <mutex>

#include "utils.h"
#include "sam_utils.h"
#include "libs/ssw_cpp.h"
#include "libs/ssw.h"

extern std::ofstream flog;
extern std::mutex log_mtx;
extern config_t config;

struct del_ins_t {
    int del, ins;

    del_ins_t() : del(0), ins(0) {};
};
del_ins_t get_dels_ins_in_first_n_chars(bam_redux_t* r, int n) {
    del_ins_t del_ins;
    int offset = 0;
    for (uint32_t c : r->cigar) {
        int len = bam_cigar_oplen(c);
        char op = bam_cigar_opchr(c);

        // since we are "unrolling" soft-clipped bases, they must be accounted for
        bool consumes_ref = (bam_cigar_type(c) & 2) || op == 'S';
        if (consumes_ref && offset + len > n) {
            len = n-offset;
        }

        if (op == 'D') {
            del_ins.del += len;
        } else if (op == 'I') {
            del_ins.ins += len;
        }

        if (consumes_ref) {
            offset += len;
            if (offset == n) break;
        }
    }
    return del_ins;
}

struct base_score_t {
    int freq = 0, qual = 0;
    char base;

    base_score_t(char base) : base(base) {}
};
bool operator < (const base_score_t& bs1, const base_score_t& bs2) {
    if (bs1.freq != bs2.freq) return bs1.freq < bs2.freq;
    return bs1.qual < bs2.qual;
}

hts_pos_t get_start_offset(bam_redux_t* r1, bam_redux_t* r2) {
    hts_pos_t offset = r2->unclipped_start() - r1->unclipped_start();
    /* suppose the difference in starting position between R1 and R2 is N, we may be tempted to align the start
     * of R2 to position N+1 of R1
     * However, if R1 has insertions or deletions in the first N bps, the difference in starting positions
     * is no longer accurate to decide where the start of R2 aligns on R1.
     * If R1 has I bps inserted in the first N bps, then R2 will align to position N+1+I.
     * Conversely, if D bps are deleted, R2 will align to position N+1-D
     */
    del_ins_t del_ins = get_dels_ins_in_first_n_chars(r1, offset);
    return offset + del_ins.ins - del_ins.del;
}

// mismatch_score, gap_open_score and gap_extend_score must be negative
int compute_read_score(bam_redux_t* r, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
	std::vector<uint32_t>& cigar = r->cigar;
	int matches = 0;
	int mismatches = r->nm;
	int score = 0;
	for (int i = 0; i < cigar.size(); i++) {
		char op = bam_cigar_opchr(cigar[i]);
		int len = bam_cigar_oplen(cigar[i]);
		if (op == 'M') matches += bam_cigar_oplen(cigar[i]);
		else if (op == 'I' || op == 'D') {
			score += gap_open_score + (len-1)*gap_extend_score;
			mismatches -= len;
		}
	}
	score += match_score*matches + mismatch_score*mismatches - match_score*mismatches;
	return score;
}

std::string build_full_consensus_seq(std::vector<bam_redux_t*>& clipped) {
    const int MAX_CONSENSUS_LEN = 100000;
    char consensus[MAX_CONSENSUS_LEN];

    std::vector<hts_pos_t> read_start_offsets, read_end_offsets;
    hts_pos_t consensus_len = 0;
    for (bam_redux_t* r : clipped) {
        hts_pos_t start_offset = get_start_offset(clipped[0], r);
        read_start_offsets.push_back(start_offset);

        hts_pos_t end_offset = start_offset + r->seq_len() - 1;
        read_end_offsets.push_back(end_offset);
        consensus_len = std::max(consensus_len, end_offset+1);
    }

    int s = 0;
    for (int i = 0; i < consensus_len; i++) {
        while (s < clipped.size() && read_end_offsets[s] < i) s++;

        base_score_t base_scores[4] = { base_score_t('A'), base_score_t('C'), base_score_t('G'), base_score_t('T') };
        for (int j = s; j < clipped.size() && read_start_offsets[j] <= i; j++) {
            if (read_end_offsets[j] < i) continue;

            char nucl = get_base(clipped[j]->seq.data(), i - read_start_offsets[j]);
            uint8_t qual = clipped[j]->qual[i - read_start_offsets[j]];
            if (nucl == 'A') {
                base_scores[0].freq++;
                if (base_scores[0].qual < qual) base_scores[0].qual = qual;
            } else if (nucl == 'C') {
                base_scores[1].freq++;
                if (base_scores[1].qual < qual) base_scores[1].qual = qual;
            } else if (nucl == 'G') {
                base_scores[2].freq++;
                if (base_scores[2].qual < qual) base_scores[2].qual = qual;
            } else if (nucl == 'T') {
                base_scores[3].freq++;
                if (base_scores[3].qual < qual) base_scores[3].qual = qual;
            }
        }

        base_score_t best_base_score = max(base_scores[0], base_scores[1], base_scores[2], base_scores[3]);

        consensus[i] = best_base_score.base;
    }
    consensus[consensus_len] = '\0';

    return std::string(consensus);
}

consensus_t* build_full_consensus(int contig_id, std::vector<bam_redux_t*>& clipped, bool left_clipped) {
    if (clipped.size() <= 2) return NULL;

    std::sort(clipped.begin(), clipped.end(), [](bam_redux_t* r1, bam_redux_t* r2) {
        return r1->unclipped_start() < r2->unclipped_start();
    });
    if (clipped[0]->unclipped_start() < 0) return NULL;

    std::string consensus_seq = build_full_consensus_seq(clipped);
    if (consensus_seq == "") return NULL;

    std::vector<bam_redux_t*> accepted_reads;
    std::vector<std::pair<int, int> > old_and_new_scores;
    for (bam_redux_t* r : clipped) {
    	int mm = 0;
        int mm_clip = 0; // mismatches in the clipped portion only

        // filter reads with too many differences from the consensus_seq
        hts_pos_t offset = get_start_offset(clipped[0], r);
        hts_pos_t clip_start = left_clipped ? 0 : r->seq_len() - r->right_clip_size;
        hts_pos_t clip_end = left_clipped ? r->left_clip_size : r->seq_len();
        for (int i = 0; i < r->seq_len(); i++) {
            if (i+offset >= consensus_seq.length()) {
                std::cerr << "WARNING: consensus_seq out of boundary." << std::endl;
            }
            if (consensus_seq[i + offset] != get_base(r->seq.data(), i)) {
                mm++;
                if (i >= clip_start && i < clip_end) mm_clip++;
            }
        }
        if (mm <= std::ceil(config.max_seq_error * consensus_seq.length()) &&
            mm_clip <= (clip_end-clip_start)/2) {
            /* this is meant to fix a corner case where the clip is very short, and completely different from
             * the consensus_seq. However, the rest of the reads is the same as the consensus_seq, therefore the mismatches
             * accumulated in the clip are not enough to discard the read, despite the read not belonging to the cluster.
             * We are being very permissive (i.e. we allow up to 50% mismatches in the clip) because
             * 1. Illumina is known to accumulate sequencing errors at the 3' tail, which is the clipped portion in about
             * 50% of the cases
             * 2. Especially for short clips, if we used config.max_seq_error as for the rest of the consensus_seq, 2-3 mismatches
             * would already discard them
             */

        	// the read should be much better when mapped to the consensus than when mapped to the reference
        	int orig_score = compute_read_score(r, 1, -4, -6, -1);
        	int new_score = (r->seq_len()-mm)*1 - mm*4;

        	if (new_score - orig_score >= config.min_score_diff) {
				old_and_new_scores.push_back({orig_score, new_score});
				accepted_reads.push_back(r);
        	}
        }
    }

    if (accepted_reads.size() <= 2) return NULL;

    hts_pos_t breakpoint = left_clipped ? INT32_MAX : 0; // the current HTS_POS_MAX does not compile on some compilers
    hts_pos_t remap_boundary = left_clipped ? consensus_t::LOWER_BOUNDARY_NON_CALCULATED : consensus_t::UPPER_BOUNDARY_NON_CALCULATED;
    int supp_clipped_reads = 0;
    uint8_t max_mapq = 0;
    std::string read_qnames;
    int i = 0;
    for (bam_redux_t* r : accepted_reads) {
        breakpoint = left_clipped ? std::min(breakpoint, r->start-1) : std::max(breakpoint, r->end-1);

        if (r->is_inter_chr() || r->is_rev() == r->is_mrev()) continue;
        // Note that if the paired reads overlap, sometimes they are both clipped at the same position -
        // i.e., they are on the same side of the deletion. We need to exclude them from remap boundary calculation
        if (left_clipped && r->is_rev() && !r->mate_left_clipped()) {
            remap_boundary = std::max(remap_boundary, r->mstart);
        } else if (!left_clipped && !r->is_rev() && !r->mate_right_clipped()) {
            remap_boundary = std::min(remap_boundary, r->start+r->isize);
        }

        if ((left_clipped && r->right_clip_size < config.min_clip_len)
         || (!left_clipped && r->left_clip_size < config.min_clip_len)) {
			supp_clipped_reads++;
			max_mapq = std::max(max_mapq, r->mapq);
            read_qnames += r->qname + "," + std::to_string(r->is_rev()) + "," + std::to_string(old_and_new_scores[i].first) + "," + std::to_string(old_and_new_scores[i].second) + ",";
        }
        i++;
    }
    if (supp_clipped_reads < 3) return NULL;

    consensus_seq = build_full_consensus_seq(accepted_reads);

    // we trim the beginning and the ending that are supported by less than 3 reads
    // we need the third smallest start offset and third highest end offset in order to trim
    // TODO: now that we break ties based on base qual, can we trust two reads?
    std::vector<hts_pos_t> start_offsets, end_offsets;
    for (bam_redux_t* r : accepted_reads) {
		hts_pos_t start_offset = get_start_offset(accepted_reads[0], r);
		start_offsets.push_back(start_offset);

		hts_pos_t end_offset = start_offset + r->seq_len();
		end_offsets.push_back(end_offset);
	}
    std::sort(start_offsets.begin(), start_offsets.end());
    std::sort(end_offsets.begin(), end_offsets.end(), std::greater<hts_pos_t>());
    hts_pos_t remove_from_start = start_offsets[1];
    hts_pos_t remove_from_end = consensus_seq.length() - end_offsets[1];

    // calculate the lenght of the clip in the consensus
    int clip_len = 0, lowq_clip_portion;
    if (left_clipped) {
    	clip_len = accepted_reads[0]->left_clip_size;
    	lowq_clip_portion = remove_from_start;
    } else {
		sort(accepted_reads.begin(), accepted_reads.end(), [](bam_redux_t* r1, bam_redux_t* r2) {
			return r1->unclipped_end() > r2->unclipped_end();
		});
    	clip_len = accepted_reads[0]->right_clip_size;
    	lowq_clip_portion = remove_from_end;
    }

    hts_pos_t start = INT32_MAX, end = 0;
    for (bam_redux_t* r : accepted_reads) {
    	start = std::min(start, r->start);
    	end = std::max(start, r->end);
    }

    consensus_t* consensus = new consensus_t(left_clipped, contig_id, start, end, breakpoint, clip_len, lowq_clip_portion, consensus_seq,
    		supp_clipped_reads, max_mapq, remap_boundary);
    log_mtx.lock();
    flog << consensus->name() << " " << read_qnames << std::endl;
    log_mtx.unlock();
    return consensus;
}

#endif /* CONSENSUS_H_ */
