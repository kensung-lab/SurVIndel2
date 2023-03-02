#ifndef SURVINDEL2_STAT_TESTS_H
#define SURVINDEL2_STAT_TESTS_H

#include <set>
#include <cmath>
#include <unordered_set>
#include <random>

#include "htslib/sam.h"
#include "sam_utils.h"
#include "ks-test.h"
#include "libs/IntervalTree.h"

std::mutex mtx2;

struct read_w_cached_info_t {
    bam1_t* read;
    hts_pos_t start, end, isize;
    int64_t as;
    int aln_len;
    bool left_clipped, right_clipped, is_rev, is_mrev;
    int references = 0;

    read_w_cached_info_t(bam1_t* read) : read(bam_dup1(read)), as(get_AS_tag(read)), start(read->core.pos), end(bam_endpos(read)),
                                         aln_len(get_aligned_portion_len(read)), left_clipped(is_left_clipped(read, 0)),
                                         right_clipped(is_right_clipped(read, 0)), is_rev(bam_is_rev(read)), is_mrev(bam_is_mrev(read)),
                                         isize(read->core.isize) {}

    ~read_w_cached_info_t() { bam_destroy1(read); }
};


int cleanup_bin_size(int read_len) {
    return 200;
}

std::vector<read_w_cached_info_t*> cleanup_reads(std::vector<read_w_cached_info_t*>& reads_w_as, hts_pos_t region_start,
                                                hts_pos_t region_end, hts_pos_t original_dup_left_bp, hts_pos_t original_dup_right_bp,
                                                stats_t stats, config_t config, std::vector<double>& max_allowed_frac_normalized,
												std::default_random_engine& rng) {
    int BIN_SIZE = cleanup_bin_size(config.read_len);
    hts_pos_t region_len = region_end - region_start;
    int n_bins = (region_len-1)/BIN_SIZE + 1;
    double desired_depth = stats.min_depth;

    // compute coverage (in term of bases) that we want to achieve for each bin
//    std::vector<uint32_t> bin_desired_cov_bases(n_bins);
//    for (int i = 0; i < n_bins-1; i++) {
//        bin_desired_cov_bases[i] = desired_depth*BIN_SIZE;
//    }
//    bin_desired_cov_bases[n_bins-1] = desired_depth*(region_len%BIN_SIZE == 0 ? BIN_SIZE : region_len%BIN_SIZE);

    std::shuffle(reads_w_as.begin(), reads_w_as.end(), rng);
    std::sort(reads_w_as.begin(), reads_w_as.end(), [](read_w_cached_info_t* r1, read_w_cached_info_t* r2) {
        return r1->aln_len-r1->as < r2->aln_len-r2->as;
    });

    auto choose_best_bin = [&region_len, &n_bins](hts_pos_t rs_in_region, hts_pos_t re_in_region, int bin_size) {
        int rs_bin = rs_in_region/bin_size, re_bin = re_in_region/bin_size;
        int best_overlap = 0, best_bin;
        for (int i = rs_bin; i <= re_bin && i < n_bins; i++) {
            int bin_overlap = overlap(rs_in_region, re_in_region, i*bin_size,
                                      std::min(hts_pos_t(i+1)*bin_size, region_len));
            if (bin_overlap > best_overlap) {
                best_overlap = bin_overlap;
                best_bin = i;
            }
        }
        return best_bin;
    };

    std::vector<read_w_cached_info_t*> solid_reads_w_as;
    std::vector<int> region_coverage(region_len+1);
    for (read_w_cached_info_t* read_w_as : reads_w_as) {
        hts_pos_t rs = read_w_as->start, re = read_w_as->end;
        if (read_w_as->left_clipped && abs(rs - original_dup_left_bp) >= 5) continue;
        if (read_w_as->right_clipped && abs(re - original_dup_right_bp) >= 5) continue;

        hts_pos_t rs_in_region = std::max(hts_pos_t(0), rs-region_start), re_in_region = std::min(re-region_start, region_len);
        if (re_in_region <= 0 || rs_in_region >= region_len) continue;

        bool is_solid = false;
        for (int i = rs_in_region; i <= re_in_region; i++) {
            if (region_coverage[i] < desired_depth) {
                is_solid = true;
                break;
            }
        }

        if (is_solid) {
            for (int i = rs_in_region; i <= re_in_region; i++) {
                region_coverage[i]++;
            }
            solid_reads_w_as.push_back(read_w_as);
        }
    }

    std::vector<std::vector<int64_t>> solid_as_diff_dists(n_bins);
    std::vector<std::vector<read_w_cached_info_t*> > binned_reads(n_bins);
    for (read_w_cached_info_t* read_w_as : solid_reads_w_as) {
        hts_pos_t rs = read_w_as->start, re = read_w_as->end;
        hts_pos_t rs_in_region = std::max(hts_pos_t(0), rs-region_start), re_in_region = std::min(re-region_start, region_len);
        int best_bin = choose_best_bin(rs_in_region, re_in_region, BIN_SIZE);
        solid_as_diff_dists[best_bin].push_back(read_w_as->aln_len-read_w_as->as);
    }
    for (read_w_cached_info_t* read_w_as : reads_w_as) {
        hts_pos_t rs = read_w_as->start, re = read_w_as->end;
        hts_pos_t rs_in_region = std::max(hts_pos_t(0), rs-region_start), re_in_region = std::min(re-region_start, region_len);
        if (re_in_region <= 0 || rs_in_region >= region_len) continue;
        int best_bin = choose_best_bin(rs_in_region, re_in_region, BIN_SIZE);
        binned_reads[best_bin].push_back(read_w_as);
    }

    std::vector<read_w_cached_info_t*> kept_reads;
    for (int i = 0; i < n_bins; i++) {
        std::vector<int64_t>& bin_as_dists = solid_as_diff_dists[i];
        int n = bin_as_dists.size();
        if (n < 3) continue;

        std::sort(bin_as_dists.begin(), bin_as_dists.end());
        int64_t median_as_diff = bin_as_dists[n/2];

        int median_or_better = 0;
        for (read_w_cached_info_t* read_w_as : binned_reads[i]) {
            if (read_w_as->aln_len-read_w_as->as <= median_as_diff) {
                median_or_better++;
            }
        }

        std::vector<int> as_diff_count(max_allowed_frac_normalized.size());
        for (read_w_cached_info_t* read_w_as : binned_reads[i]) {
            if (read_w_as->left_clipped && abs(read_w_as->start - original_dup_left_bp) >= 5) continue;
            if (read_w_as->right_clipped && abs(read_w_as->end - original_dup_right_bp) >= 5) continue;

            int64_t as_diff = read_w_as->aln_len - read_w_as->as;
            int64_t shifted_as_diff = std::max(int64_t(0), as_diff-median_as_diff);
            if (shifted_as_diff >= max_allowed_frac_normalized.size()) continue;
            if (as_diff_count[shifted_as_diff]+1 > max_allowed_frac_normalized[shifted_as_diff] * median_or_better) continue;

            kept_reads.push_back(read_w_as);
            as_diff_count[shifted_as_diff]++;
        }
    }

    return kept_reads;
}

void compute_cleaned_up_depth(duplication_t* dup, open_samFile_t* bam_file,
       std::vector<read_w_cached_info_t*>& lf_reads, std::vector<read_w_cached_info_t*>& ldup_reads,
       std::vector<read_w_cached_info_t*>& rdup_reads, std::vector<read_w_cached_info_t*>& rf_reads,
       stats_t& stats, config_t& config, std::vector<double>& max_allowed_frac_normalized, std::string workdir) {

	const int FLANKING_SIZE = 5000, INDEL_TESTED_REGION_SIZE = 10000;

    hts_pos_t left_flanking_start = dup->start-FLANKING_SIZE, left_flanking_end = dup->start;
    hts_pos_t dup_left_start = dup->start, dup_left_end = std::min(dup->end, dup->start+INDEL_TESTED_REGION_SIZE);
    hts_pos_t dup_right_start = std::max(dup->start, dup->end-INDEL_TESTED_REGION_SIZE), dup_right_end = dup->end;
    hts_pos_t right_flanking_start = dup->end, right_flanking_end = dup->end+FLANKING_SIZE;

    // cleanup reads
//    std::vector<read_w_cached_info_t*> left_flanking_reads =
//            cleanup_reads(lf_reads, left_flanking_start, left_flanking_end, dup->original_start, dup->original_end,
//                          stats, config, max_allowed_frac_normalized);
//    std::vector<read_w_cached_info_t*> dup_left_reads =
//            cleanup_reads(ldup_reads, dup_left_start, dup_left_end, dup->original_start, dup->original_end, stats, config,
//                          max_allowed_frac_normalized);
//    std::vector<read_w_cached_info_t*> dup_right_reads =
//            cleanup_reads(rdup_reads, dup_right_start, dup_right_end, dup->original_start, dup->original_end, stats, config,
//                          max_allowed_frac_normalized);
//    std::vector<read_w_cached_info_t*> right_flanking_reads =
//            cleanup_reads(rf_reads, right_flanking_start, right_flanking_end, dup->original_start, dup->original_end,
//                          stats, config,max_allowed_frac_normalized);

    std::set<read_w_cached_info_t*> kept_reads;
    kept_reads.insert(lf_reads.begin(), lf_reads.end());
    kept_reads.insert(ldup_reads.begin(), ldup_reads.end());
    kept_reads.insert(rdup_reads.begin(), rdup_reads.end());
    kept_reads.insert(rf_reads.begin(), rf_reads.end());

    uint32_t left_flanking_coverage[FLANKING_SIZE], dup_left_coverage[INDEL_TESTED_REGION_SIZE];
    uint32_t dup_right_coverage[INDEL_TESTED_REGION_SIZE], right_flanking_coverage[FLANKING_SIZE];
    memset(left_flanking_coverage, 0, FLANKING_SIZE * sizeof(left_flanking_coverage[0]));
    memset(dup_left_coverage, 0, INDEL_TESTED_REGION_SIZE * sizeof(dup_left_coverage[0]));
    memset(dup_right_coverage, 0, INDEL_TESTED_REGION_SIZE * sizeof(dup_right_coverage[0]));
    memset(right_flanking_coverage, 0, FLANKING_SIZE * sizeof(right_flanking_coverage[0]));
    for (read_w_cached_info_t* read_w_as : kept_reads) {
        hts_pos_t rs = read_w_as->start, re = read_w_as->end;

        // note: this may overestimate the coverage, for example if the read has deletions. But it is much simpler and faster than considering cigar
        hts_pos_t b = std::max(hts_pos_t(0), rs-left_flanking_start), e = std::min(hts_pos_t(FLANKING_SIZE), re-left_flanking_start);
        for (hts_pos_t i = b; i < e; i++) {
            left_flanking_coverage[i]++;
        }

        b = std::max(hts_pos_t(0), rs-dup_left_start), e = std::min(hts_pos_t(INDEL_TESTED_REGION_SIZE), re-dup_left_start);
        for (hts_pos_t i = b; i < e; i++) {
            dup_left_coverage[i]++;
        }

        b = std::max(hts_pos_t(0), rs-dup_right_start), e = std::min(hts_pos_t(INDEL_TESTED_REGION_SIZE), re-dup_right_start);
        for (hts_pos_t i = b; i < e; i++) {
            dup_right_coverage[i]++;
        }

        b = std::max(hts_pos_t(0), rs-right_flanking_start), e = std::min(hts_pos_t(FLANKING_SIZE), re-right_flanking_start);
        for (hts_pos_t i = b; i < e; i++) {
            right_flanking_coverage[i]++;
        }
    }

    // trying median
    hts_pos_t dup_tested_len = std::min(dup->end-dup->start, hts_pos_t(INDEL_TESTED_REGION_SIZE));
    std::sort(left_flanking_coverage, left_flanking_coverage+FLANKING_SIZE, std::greater<uint32_t>());
    int end = std::find(left_flanking_coverage, left_flanking_coverage+FLANKING_SIZE, 0)-left_flanking_coverage;
    dup->med_left_flanking_cov = dup->left_flanking_cov = left_flanking_coverage[end/2];
    std::sort(dup_left_coverage, dup_left_coverage+dup_tested_len, std::greater<uint32_t>());
    dup->med_indel_left_cov = dup->indel_left_cov = dup_left_coverage[dup_tested_len/2];
    std::sort(dup_right_coverage, dup_right_coverage+dup_tested_len, std::greater<uint32_t>());
    dup->med_indel_right_cov = dup->indel_right_cov = dup_right_coverage[dup_tested_len/2];
    std::sort(right_flanking_coverage, right_flanking_coverage+FLANKING_SIZE, std::greater<uint32_t>());
    end = std::find(right_flanking_coverage, right_flanking_coverage+FLANKING_SIZE, 0)-right_flanking_coverage;
    dup->med_right_flanking_cov = dup->right_flanking_cov = right_flanking_coverage[end/2];

    int ow_pairs = 0;

    for (read_w_cached_info_t* r : ldup_reads) {
        hts_pos_t mate_end = r->start + r->isize;
        if (r->start >= dup->original_start && r->end-dup->original_start <= config.max_is && r->is_rev && !r->is_mrev
        && mate_end <= dup->original_end && dup->original_end-mate_end <= config.max_is && r->isize > 0) {
            ow_pairs++;
        }
    }
    dup->disc_pairs = ow_pairs;
}


void clear_rcis(std::vector<read_w_cached_info_t*>& testable_dups_lf_reads, std::vector<read_w_cached_info_t*>& testable_dups_ldup_reads,
                std::vector<read_w_cached_info_t*>& testable_dups_rdup_reads, std::vector<read_w_cached_info_t*>& testable_dups_rf_reads) {
    for (read_w_cached_info_t* r : testable_dups_lf_reads) {
        r->references--;
        if (r->references == 0) delete r;
    }
    for (read_w_cached_info_t* r : testable_dups_ldup_reads) {
        r->references--;
        if (r->references == 0) delete r;
    }
    for (read_w_cached_info_t* r : testable_dups_rdup_reads) {
        r->references--;
        if (r->references == 0) delete r;
    }
    for (read_w_cached_info_t* r : testable_dups_rf_reads) {
        r->references--;
        if (r->references == 0) delete r;
    }
    std::vector<read_w_cached_info_t*>().swap(testable_dups_lf_reads);
    std::vector<read_w_cached_info_t*>().swap(testable_dups_ldup_reads);
    std::vector<read_w_cached_info_t*>().swap(testable_dups_rdup_reads);
    std::vector<read_w_cached_info_t*>().swap(testable_dups_rf_reads);
}

void depth_filter_dup_w_cleanup(std::string contig_name, std::vector<duplication_t*>& duplications,
                                open_samFile_t* bam_file, stats_t stats, config_t config,
                                std::vector<double>& max_allowed_frac_normalized, std::string workdir) {

	if (duplications.empty()) return;

	const int FLANKING_SIZE = 5000, INDEL_TESTED_REGION_SIZE = 10000;

	std::default_random_engine rng(config.seed);

	std::vector<hts_pair_pos_t> cleanup_ivals;
	for (duplication_t* dup : duplications) {
		cleanup_ivals.push_back({std::max(dup->start-FLANKING_SIZE, hts_pos_t(0)), std::min(dup->start+INDEL_TESTED_REGION_SIZE, dup->end)});
		cleanup_ivals.push_back({std::max(dup->end-INDEL_TESTED_REGION_SIZE, dup->start), dup->end+FLANKING_SIZE});
	}
	std::sort(cleanup_ivals.begin(), cleanup_ivals.end(), [](hts_pair_pos_t p1, hts_pair_pos_t p2) {
		return p1.beg < p2.beg;
	});

	std::vector<hts_pair_pos_t> merged_cleanup_ivals;
	hts_pos_t curr_start = cleanup_ivals[0].beg, curr_end = cleanup_ivals[0].end;
	for (hts_pair_pos_t cleanup_ival : cleanup_ivals) {
		if (cleanup_ival.beg > curr_end) {
			merged_cleanup_ivals.push_back({curr_start, curr_end});
			curr_start = cleanup_ival.beg, curr_end = cleanup_ival.end;
		}
		curr_end = std::max(curr_end, cleanup_ival.end);
	}
	merged_cleanup_ivals.push_back({curr_start, curr_end});

	bam1_t* read = bam_init1();
	std::vector<std::vector<read_w_cached_info_t*> > clean_reads;
	for (hts_pair_pos_t merged_cleanup_ival : merged_cleanup_ivals) {
		std::vector<read_w_cached_info_t*> reads;

		std::stringstream ss;
		ss << contig_name << ":" << merged_cleanup_ival.beg << "-" << merged_cleanup_ival.end;
		hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, ss.str().c_str());
		while (sam_itr_next(bam_file->file, iter, read) >= 0) {
			if (is_unmapped(read) || is_mate_unmapped(read)) continue;
			if (read->core.tid != read->core.mtid || read->core.qual < 20) continue; // do not count low qual reads and inter-chr pairs

			reads.push_back(new read_w_cached_info_t(read));
		}

		std::vector<read_w_cached_info_t*> curr_clean_reads = cleanup_reads(reads, merged_cleanup_ival.beg, merged_cleanup_ival.end,
				merged_cleanup_ival.beg, merged_cleanup_ival.end, stats, config, max_allowed_frac_normalized, rng);
		std::vector<read_w_cached_info_t*> curr_clean_reads_copy;
		for (read_w_cached_info_t* rci : curr_clean_reads) {
			curr_clean_reads_copy.push_back(new read_w_cached_info_t(rci->read));
		}
		clean_reads.push_back(curr_clean_reads_copy);

		for (read_w_cached_info_t* rci : reads) delete rci;
	}

	std::vector<Interval<int> > ivals;
	for (int i = 0; i < merged_cleanup_ivals.size(); i++) {
		ivals.push_back(Interval<int>(merged_cleanup_ivals[i].beg, merged_cleanup_ivals[i].end, i));
	}
	IntervalTree<int> ival_tree(ivals);

	for (duplication_t* dup : duplications) {
		std::vector<Interval<int> > c_ivals = ival_tree.findOverlapping(std::max(dup->start-FLANKING_SIZE, hts_pos_t(0)),
				std::min(dup->start+INDEL_TESTED_REGION_SIZE, dup->end));
		int idx = c_ivals[0].value;

		std::vector<read_w_cached_info_t*> dups_lf_reads, dups_ldup_reads;
		std::vector<read_w_cached_info_t*> dups_rdup_reads, dups_rf_reads;
		hts_pos_t left_flanking_start = dup->start - FLANKING_SIZE, left_flanking_end = dup->start;
		hts_pos_t dup_left_start = dup->start, dup_left_end = std::min(dup->end, dup->start+INDEL_TESTED_REGION_SIZE);
		hts_pos_t dup_right_start = std::max(dup->start, dup->end-INDEL_TESTED_REGION_SIZE), dup_right_end = dup->end;
		hts_pos_t right_flanking_start = dup->end, right_flanking_end = dup->end + FLANKING_SIZE;

		for (read_w_cached_info_t* rci : clean_reads[idx]) {
			bam1_t* read = rci->read;
			hts_pos_t rs = read->core.pos, re = bam_endpos(read);
			bool overlaps_lf = overlap(rs, re, left_flanking_start, left_flanking_end) > 0;
			bool overlaps_ldup = overlap(rs, re, dup_left_start, dup_left_end) > 0;
			bool overlaps_rdup = overlap(rs, re, dup_right_start, dup_right_end) > 0;
			bool overlaps_rf = overlap(rs, re, right_flanking_start, right_flanking_end) > 0;
			if (overlaps_lf || overlaps_ldup || overlaps_rdup || overlaps_rf) {
				if (overlaps_lf) {
					dups_lf_reads.push_back(rci);
				}
				if (overlaps_ldup) {
					dups_ldup_reads.push_back(rci);
				}
				if (overlaps_rdup) {
					dups_rdup_reads.push_back(rci);
				}
				if (overlaps_rf) {
					dups_rf_reads.push_back(rci);
				}
			}
		}

		compute_cleaned_up_depth(dup, bam_file, dups_lf_reads, dups_ldup_reads, dups_rdup_reads, dups_rf_reads, stats,
				config, max_allowed_frac_normalized, workdir);
	}

	for (auto& v : clean_reads) {
		for (read_w_cached_info_t* rci : v) {
			delete rci;
		}
	}

    // we have our population-wise distribution of AS differences
    // any read with AS >= the median for its best bin is considered AS diff = 0
    // from there, we compute how many reads we allow for each AS diff > 0, based on the population-wise distribution
//    std::vector<char*> regions;
//    for (duplication_t* dup : duplications) {
//        std::stringstream ss;
//        ss << contig_name << ":" << dup->start-FLANKING_SIZE << "-" << std::min(dup->start+INDEL_TESTED_REGION_SIZE, dup->end);
//        char* region = new char[1000];
//        strcpy(region, ss.str().c_str());
//        regions.push_back(region);
//
//        ss.str(std::string());
//        ss << contig_name << ":" << std::max(dup->end-INDEL_TESTED_REGION_SIZE, dup->start) << "-" << dup->end+FLANKING_SIZE;
//        region = new char[1000];
//        strcpy(region, ss.str().c_str());
//        regions.push_back(region);
//    }
//    std::sort(duplications.begin(), duplications.end(), [](const duplication_t* dup1, const duplication_t* dup2) {
//        return dup1->start < dup2->start;
//    });
//
//    std::vector<std::pair<duplication_t*, int> > duplications_endsorted(duplications.size());
//    for (int i = 0; i < duplications.size(); i++) {
//        duplications_endsorted[i] = std::make_pair(duplications[i], i);
//    }
//    std::sort(duplications_endsorted.begin(), duplications_endsorted.end(),
//              [] (std::pair<duplication_t*, int>& p1, std::pair<duplication_t*, int>& p2)
//              { return p1.first->end < p2.first->end; });
//
//
//    std::vector<std::vector<read_w_cached_info_t*> > dups_lf_reads(duplications.size()), dups_ldup_reads(duplications.size());
//    std::vector<std::vector<read_w_cached_info_t*> > dups_rdup_reads(duplications.size()), dups_rf_reads(duplications.size());
//    curr_pos = 0;
//    int processed_dups = 0;
//
//    for (bam1_t* read : reads) {
//        hts_pos_t rs = read->core.pos, re = bam_endpos(read);
//        while (curr_pos < duplications.size() && rs > duplications[curr_pos]->end+FLANKING_SIZE) curr_pos++;
//
//        while (processed_dups < duplications.size() && rs > duplications_endsorted[processed_dups].first->end+FLANKING_SIZE) {
//            int idx = duplications_endsorted[processed_dups].second;
//            compute_cleaned_up_depth(duplications[idx], bam_file, dups_lf_reads[idx],
//                                     dups_ldup_reads[idx], dups_rdup_reads[idx],
//                                     dups_rf_reads[idx], stats, config, max_allowed_frac_normalized, workdir);
//            clear_rcis(dups_lf_reads[idx],dups_ldup_reads[idx],dups_rdup_reads[idx],dups_rf_reads[idx]);
//            processed_dups++;
//        }
//
//        read_w_cached_info_t* rci = new read_w_cached_info_t(read);
//        bool delete_rci = true;
//        for (int i = curr_pos; i < duplications.size() && re > duplications[i]->start-FLANKING_SIZE; i++) {
//            duplication_t* dup = duplications[i];
//            hts_pos_t left_flanking_start = dup->start - FLANKING_SIZE, left_flanking_end = dup->start;
//            hts_pos_t dup_left_start = dup->start, dup_left_end = std::min(dup->end, dup->start+INDEL_TESTED_REGION_SIZE);
//            hts_pos_t dup_right_start = std::max(dup->start, dup->end-INDEL_TESTED_REGION_SIZE), dup_right_end = dup->end;
//            hts_pos_t right_flanking_start = dup->end, right_flanking_end = dup->end + FLANKING_SIZE;
//
//            bool overlaps_lf = overlap(rs, re, left_flanking_start, left_flanking_end) > 0;
//            bool overlaps_ldup = overlap(rs, re, dup_left_start, dup_left_end) > 0;
//            bool overlaps_rdup = overlap(rs, re, dup_right_start, dup_right_end) > 0;
//            bool overlaps_rf = overlap(rs, re, right_flanking_start, right_flanking_end) > 0;
//
//            if (overlaps_lf || overlaps_ldup || overlaps_rdup || overlaps_rf) {
//                if (overlaps_lf) {
//                    dups_lf_reads[i].push_back(rci);
//                }
//                if (overlaps_ldup) {
//                    dups_ldup_reads[i].push_back(rci);
//                }
//                if (overlaps_rdup) {
//                    dups_rdup_reads[i].push_back(rci);
//                }
//                if (overlaps_rf) {
//                    dups_rf_reads[i].push_back(rci);
//                }
//                rci->references += overlaps_lf + overlaps_ldup + overlaps_rdup + overlaps_rf;
//                delete_rci = false;
//            }
//        }
//        if (delete_rci) delete rci;
//    }
//    while (processed_dups < duplications.size()) {
//        int idx = duplications_endsorted[processed_dups].second;
//        compute_cleaned_up_depth(duplications[idx], bam_file, dups_lf_reads[idx],
//                                 dups_ldup_reads[idx], dups_rdup_reads[idx],
//                dups_rf_reads[idx], stats, config, max_allowed_frac_normalized, workdir);
//        clear_rcis(dups_lf_reads[idx],dups_ldup_reads[idx],dups_rdup_reads[idx],dups_rf_reads[idx]);
//        processed_dups++;
//    }
//
//    for (char* region : regions) {
//        delete[] region;
//    }
//    hts_itr_destroy(iter);
//    bam_destroy1(read);
}

void depth_filter_indel(std::string contig_name, std::vector<indel_t*>& indels, open_samFile_t* bam_file, config_t& config) {

	const int FLANKING_SIZE = 5000, INDEL_TESTED_REGION_SIZE = 10000;

	std::sort(indels.begin(), indels.end(), [](const indel_t* i1, const indel_t* i2) {
		return i1->start < i2->start;
	});

    std::vector<char*> regions;
    for (indel_t* indel : indels) {
        std::stringstream ss;
        ss << contig_name << ":" << std::max(hts_pos_t(0), indel->start - FLANKING_SIZE) << "-" << std::min(indel->start + INDEL_TESTED_REGION_SIZE, indel->end);
        char* region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);

        ss.str(std::string());
        ss << contig_name << ":" << std::max(indel->end - INDEL_TESTED_REGION_SIZE, indel->start) << "-" << indel->end + FLANKING_SIZE;
        region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);
    }
    std::sort(indels.begin(), indels.end(), [](const indel_t* indel1, const indel_t* indel2) {
        return indel1->start < indel2->start;
    });

    std::vector<int64_t> flanking_left_cov(indels.size()), del_left_cov(indels.size());
    std::vector<int64_t> del_right_cov(indels.size()), flanking_right_cov(indels.size());

    std::vector<std::vector<int64_t> > _flanking_left_cov(indels.size(), std::vector<int64_t>(FLANKING_SIZE)), _indel_left_cov(indels.size(), std::vector<int64_t>(INDEL_TESTED_REGION_SIZE));
	std::vector<std::vector<int64_t> > _indel_right_cov(indels.size(), std::vector<int64_t>(INDEL_TESTED_REGION_SIZE)), _flanking_right_cov(indels.size(), std::vector<int64_t>(FLANKING_SIZE));
	std::vector<std::vector<int64_t> > _left_cluster_cov(indels.size(), std::vector<int64_t>(config.max_is+1)), _right_cluster_cov(indels.size(), std::vector<int64_t>(config.max_is+1));

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    int curr_pos = 0;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
        if (!is_samechr(read)) continue; // do not count low qual reads and inter-chr pairs

        hts_pos_t rs = read->core.pos, re = bam_endpos(read);
        while (curr_pos < indels.size() && indels[curr_pos]->end + FLANKING_SIZE < rs) curr_pos++;

        if (curr_pos == indels.size()) break;

        for (int i = curr_pos; i < indels.size() && indels[i]->start - FLANKING_SIZE < re; i++) {
            indel_t* indel = indels[i];
			hts_pos_t indel_tested_len = std::min(indel->end-indel->start, hts_pos_t(INDEL_TESTED_REGION_SIZE));

			hts_pos_t left_flanking_start = std::max(hts_pos_t(0), indel->start-FLANKING_SIZE);
			hts_pos_t left_flanking_len = indel->start - left_flanking_start;
			hts_pos_t indel_lcov_start = indel->start;
			hts_pos_t indel_rcov_start = std::max(indel->start, indel->end-INDEL_TESTED_REGION_SIZE);
			hts_pos_t right_flanking_start = indel->end;

			hts_pos_t b = std::max(hts_pos_t(0), rs-indel->rc_anchor_start), e = std::min(re-indel->rc_anchor_start, indel->start-indel->rc_anchor_start);
			for (hts_pos_t j = b; j < e; j++) {
				_left_cluster_cov[i][j]++;
			}
			b = std::max(hts_pos_t(0), rs-indel->end), e = std::min(re-indel->end, indel->lc_anchor_end-indel->end);
			for (hts_pos_t j = b; j < e; j++) {
				_right_cluster_cov[i][j]++;
			}

			if (read->core.qual < 20) continue;

			b = std::max(hts_pos_t(0), rs-left_flanking_start), e = std::min(left_flanking_len, re-left_flanking_start);
			for (hts_pos_t j = b; j < e; j++) {
				_flanking_left_cov[i][j]++;
			}

			b = std::max(hts_pos_t(0), rs-indel_lcov_start), e = std::min(indel_tested_len, re-indel_lcov_start);
			for (hts_pos_t j = b; j < e; j++) {
				_indel_left_cov[i][j]++;
			}

			b = std::max(hts_pos_t(0), rs-indel_rcov_start), e = std::min(indel_tested_len, re-indel_rcov_start);
			for (hts_pos_t j = b; j < e; j++) {
				_indel_right_cov[i][j]++;
			}

			b = std::max(hts_pos_t(0), rs-right_flanking_start), e = std::min(hts_pos_t(FLANKING_SIZE), re-right_flanking_start);
			for (hts_pos_t j = b; j < e; j++) {
				_flanking_right_cov[i][j]++;
			}

            flanking_left_cov[i] += overlap(indel->start - FLANKING_SIZE, indel->start, rs, re);
            del_left_cov[i] += overlap(indel->start, std::min(indel->start + INDEL_TESTED_REGION_SIZE, indel->end), rs, re);
            del_right_cov[i] += overlap(std::max(indel->end - INDEL_TESTED_REGION_SIZE, indel->start), indel->end, rs, re);
            flanking_right_cov[i] += overlap(indel->end, indel->end + FLANKING_SIZE, rs, re);
        }
    }

    for (int i = 0; i < indels.size(); i++) {
        indel_t* indel = indels[i];
		hts_pos_t indel_tested_len = std::min(indel->end-indel->start, hts_pos_t(INDEL_TESTED_REGION_SIZE));

		indel->left_flanking_cov = flanking_left_cov[i] / FLANKING_SIZE;
        indel->indel_left_cov = del_left_cov[i] / indel_tested_len;
        indel->indel_right_cov = del_right_cov[i] / indel_tested_len;
        indel->right_flanking_cov = flanking_right_cov[i] / FLANKING_SIZE;

        // trying median
		std::sort(_flanking_left_cov[i].begin(), _flanking_left_cov[i].end());
		indel->med_left_flanking_cov = _flanking_left_cov[i][_flanking_left_cov[i].size()/2];

		std::sort(_indel_left_cov[i].begin(), _indel_left_cov[i].end(), std::greater<int64_t>());
		indel->med_indel_left_cov = _indel_left_cov[i][indel_tested_len/2];

		std::sort(_indel_right_cov[i].begin(), _indel_right_cov[i].end(), std::greater<int64_t>());
		indel->med_indel_right_cov = _indel_right_cov[i][indel_tested_len/2];

		std::sort(_flanking_right_cov[i].begin(), _flanking_right_cov[i].end());
		indel->med_right_flanking_cov = _flanking_right_cov[i][_flanking_right_cov[i].size()/2];

		if (indel->indel_type() == "DEL") {
			_left_cluster_cov[i].resize(indel->start-indel->rc_anchor_start);
			std::sort(_left_cluster_cov[i].begin(), _left_cluster_cov[i].end());
			indel->med_left_cluster_cov = _left_cluster_cov[i][_left_cluster_cov[i].size()/2];

			_right_cluster_cov[i].resize(indel->lc_anchor_end-indel->end);
			std::sort(_right_cluster_cov[i].begin(), _right_cluster_cov[i].end());
			indel->med_right_cluster_cov = _right_cluster_cov[i][_right_cluster_cov[i].size()/2];
		}
    }

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void depth_filter_del(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file,
                      int min_size_for_depth_filtering, config_t& config) {
    std::vector<indel_t*> testable_indels(deletions.begin(), deletions.end());
    depth_filter_indel(contig_name, testable_indels, bam_file, config);
}
void depth_filter_dup(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file,
                      int min_size_for_depth_filtering, config_t& config) {
    std::vector<indel_t*> testable_indels(duplications.begin(), duplications.end());
    depth_filter_indel(contig_name, testable_indels, bam_file, config);
}

int find_smallest_range_start(std::vector<int>& v, int range_size, int& min_cum) {
	int cum = 0;
	for (int j = 0; j < range_size && j < v.size(); j++) {
		cum += v[j];
	}
	int best_start = 0, best_cum = cum;
	for (int j = range_size; j < v.size(); j++) {
		cum = cum - v[j-range_size] + v[j];
		if (cum < best_cum) {
			best_cum = cum;
			best_start = j - range_size + 1;
		}
	}
	min_cum = best_cum;
	return best_start;
}

void calculate_confidence_interval_size(std::string contig_name, std::vector<double>& global_crossing_isize_dist,
										std::vector<uint32_t>& median_crossing_count_geqi_by_isize,
										std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, config_t config, stats_t stats) {

    std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
        return (d1->start+d1->end)/2 < (d2->start+d2->end)/2;
    });

    std::vector<hts_pos_t> midpoints, sizes;
    std::vector<uint64_t> sums(deletions.size()), sq_sums(deletions.size());
    std::vector<uint32_t> ns(deletions.size());
    std::vector<char*> regions;
    std::vector<std::pair<hts_pos_t, hts_pos_t> > regions_coos;
    for (deletion_t* deletion : deletions) {
        hts_pos_t midpoint = (deletion->start+deletion->end)/2, size = deletion->end-deletion->start;
        midpoints.push_back(midpoint);
        sizes.push_back(size);

        regions_coos.push_back({deletion->start-config.max_is, deletion->start});
        regions_coos.push_back({midpoint-config.max_is, midpoint});
    }
	std::sort(regions_coos.begin(), regions_coos.end());
    for (auto coos : regions_coos) {
    	std::stringstream ss;
    	ss << contig_name << ":" << coos.first << "-" << coos.second;
    	char* region = new char[1000];
    	strcpy(region, ss.str().c_str());
		regions.push_back(region);
    }

	int curr_pos = 0;
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    std::vector<std::vector<double> > local_dists(deletions.size());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {

    	if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read) || read->core.qual < 20) continue;
        if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0) continue;

        while (curr_pos < deletions.size() && midpoints[curr_pos] < read->core.pos) curr_pos++;

        hts_pos_t start = read->core.pos + read->core.l_qseq/2;
        hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;
        for (int i = curr_pos; i < deletions.size() && midpoints[i] < read->core.pos+read->core.isize; i++) {
            if (start <= midpoints[i] && midpoints[i] <= end && read->core.isize <= config.max_is+sizes[i]) {
                sums[i] += read->core.isize;
                sq_sums[i] += read->core.isize*read->core.isize;
                ns[i]++;
                local_dists[i].push_back(read->core.isize);
            }
        }
    }

    for (int i = 0; i < deletions.size(); i++) {
        deletion_t* del = deletions[i];
        uint32_t n = ns[i];
        uint64_t sum = sums[i], sq_sum = sq_sums[i];
        if (n >= 4) {
            int avg_is = sum/n;
            int var_is = (sq_sum - sum*sum/n)/(n-1);
            int confidence_ival = 2.576 * sqrt(var_is/n);
            del->max_conf_size = avg_is - stats.pop_avg_crossing_is + confidence_ival;
            if (!global_crossing_isize_dist.empty()) {
				double p_val = ks_test(global_crossing_isize_dist, local_dists[i]);
				del->ks_pval = p_val;
				if (!del->imprecise()) continue;

				int est_size = avg_is - stats.pop_avg_crossing_is;

				// compute depth base by base in the imprecise deleted regions
				// TODO: there are some faster data structures out there for this - e.g., Fenwick tree
				hts_pos_t range_start = del->start - config.read_len, range_end = del->end + config.read_len;
				if (est_size > range_end-range_start || est_size < config.min_sv_size) continue; // estimated size is grossly off - ignore

				std::vector<int> depth_by_base(range_end-range_start+1);
				char region[1000];
				sprintf(region, "%s:%d-%d", contig_name.c_str(), (int) range_start, (int) range_end);
				iter = sam_itr_querys(bam_file->idx, bam_file->header, region);
				while (sam_itr_next(bam_file->file, iter, read) >= 0) {
					if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read) || read->core.qual < 0) continue;
					if (!is_samechr(read) || is_samestr(read) || read->core.isize < -config.max_is || read->core.isize > config.max_is) continue;

					int start = std::max(hts_pos_t(0), read->core.pos-range_start);
					int end = std::min(bam_endpos(read)-range_start, range_end-range_start);
					for (int i = start; i <= end; i++) depth_by_base[i]++;
				}

				// find window of est_size bp with minimum depth
				int min_cum_depth;
				int best_start = find_smallest_range_start(depth_by_base, est_size, min_cum_depth);
				int del_cov = min_cum_depth/est_size;

				std::vector<double> homalt_local_dist, het_local_dist;
				for (int i = 0; i < global_crossing_isize_dist.size(); i++) {
					double d = global_crossing_isize_dist[i];
					homalt_local_dist.push_back(d+est_size);
					if (i%2) het_local_dist.push_back(d+2*est_size);
					else het_local_dist.push_back(d);
				}

				int homalt_evidence = 0, het_evidence = 0;

				// whether the distribution supports a het or a homalt deletion
				double homalt_pval = ks_test(local_dists[i], homalt_local_dist);
				double het_pval = ks_test(local_dists[i], het_local_dist);
				if (homalt_pval < het_pval) { // deletion is het
					het_evidence++;
				} else { // deletion is homalt
					homalt_evidence++;
				}

				// whether the discordant pairs are closer to what we expect with a het or a hom alt deletion
				// TODO: this is better than random but not very accurate. If we could have a test on standard deviations (het will have larger
				// than the global, hom alt should be the same) as like Bartlett as the third vote it might be better
				int homalt_idx = std::max(0, config.max_is - est_size);
				int het_idx = std::max(0, config.max_is - 2*est_size);
				int dev_if_homalt = abs((int) (del->disc_pairs-median_crossing_count_geqi_by_isize[homalt_idx]));
				int dev_if_het = abs((int) (del->disc_pairs-median_crossing_count_geqi_by_isize[het_idx]/2));
				if (dev_if_het <= dev_if_homalt) {
					het_evidence++;
				} else { // deletion is homalt
					homalt_evidence++;
				}

				// depth supports hom alt or het
				if (del->med_left_flanking_cov*0.25 < del_cov || del->med_right_flanking_cov*0.25 < del_cov) {
					het_evidence++;
				} else {
					homalt_evidence++;
				}

				if (het_evidence > homalt_evidence) {
					del->genotype = "0/1";
				} else {
					del->genotype = "1/1";
				}

				// if deletion is estimated het, double the size (unless it is bigger than the original range - then assume GT is wrong)
				if (del->genotype == "0/1" && est_size*2 <= range_end-range_start) {
					est_size *= 2;
					best_start = find_smallest_range_start(depth_by_base, est_size, min_cum_depth);
				}

				int new_start = range_start + best_start, new_end = new_start + est_size;
				if (new_start >= del->start-config.read_len/2 && new_end <= del->end+config.read_len) {
					del->original_range = std::to_string(del->start) + "-" + std::to_string(del->end);
					del->start = new_start; del->end = new_end;
				}
            }
        }
    }

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void calculate_cluster_region_disc(std::string contig_name, std::vector<deletion_t*> deletions, open_samFile_t* bam_file) {

	std::vector<char*> l_cluster_regions, r_cluster_regions;
	for (deletion_t* deletion : deletions) {
		std::stringstream ss;
		ss << contig_name << ":" << deletion->rc_anchor_start << "-" << deletion->start;
		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		l_cluster_regions.push_back(region);

		ss.str("");
		ss << contig_name << ":" << deletion->end << "-" << deletion->lc_anchor_end;
		region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		r_cluster_regions.push_back(region);
	}

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->rc_anchor_start < d2->rc_anchor_start;
	});

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, l_cluster_regions.data(), l_cluster_regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < deletions.size() && deletions[curr_pos]->start < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_outward(read)) {
			for (int i = curr_pos; i < deletions.size() && deletions[i]->rc_anchor_start <= bam_endpos(read); i++) {
				if (read->core.pos <= deletions[i]->start) deletions[i]->l_cluster_region_disc_pairs++;
			}
		}
	}

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->end < d2->end;
	});

	curr_pos = 0;
	iter = sam_itr_regarray(bam_file->idx, bam_file->header, r_cluster_regions.data(), r_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < deletions.size() && deletions[curr_pos]->lc_anchor_end < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_outward(read)) {
			for (int i = curr_pos; i < deletions.size() && deletions[i]->end <= bam_endpos(read); i++) {
				if (read->core.pos <= deletions[i]->lc_anchor_end) deletions[i]->r_cluster_region_disc_pairs++;
			}
		}
	}

	for (char* region : l_cluster_regions) {
		delete[] region;
	}
	for (char* region : r_cluster_regions) {
		delete[] region;
	}
}

void calculate_cluster_region_disc(std::string contig_name, std::vector<duplication_t*> duplications, open_samFile_t* bam_file, config_t& config) {

	std::vector<char*> lc_cluster_regions, rc_cluster_regions;
	for (duplication_t* duplication : duplications) {
		std::stringstream ss;
		ss << contig_name << ":" << duplication->rc_anchor_start << "-" << duplication->end;
		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		rc_cluster_regions.push_back(region);

		ss.str("");
		ss << contig_name << ":" << duplication->start << "-" << duplication->lc_anchor_end;
		region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		lc_cluster_regions.push_back(region);
	}

	std::sort(duplications.begin(), duplications.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->rc_anchor_start < d2->rc_anchor_start;
	});

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, rc_cluster_regions.data(), rc_cluster_regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < duplications.size() && duplications[curr_pos]->end < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_long(read, config.max_is)) {
			for (int i = curr_pos; i < duplications.size() && duplications[i]->rc_anchor_start <= bam_endpos(read); i++) {
				if (read->core.pos <= duplications[i]->end) {
					duplications[i]->rc_cluster_region_disc_pairs++;
				}
			}
		}
	}

	std::sort(duplications.begin(), duplications.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->start < d2->start;
	});

	curr_pos = 0;
	iter = sam_itr_regarray(bam_file->idx, bam_file->header, rc_cluster_regions.data(), rc_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < duplications.size() && duplications[curr_pos]->lc_anchor_end < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_long(read, config.max_is)) {
			for (int i = curr_pos; i < duplications.size() && duplications[i]->start <= bam_endpos(read); i++) {
				if (read->core.pos <= duplications[i]->lc_anchor_end) duplications[i]->lc_cluster_region_disc_pairs++;
			}
		}
	}

	for (char* region : lc_cluster_regions) {
		delete[] region;
	}
	for (char* region : rc_cluster_regions) {
		delete[] region;
	}
}

void calculate_ptn_ratio(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, config_t config) {

	if (deletions.empty()) return;

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return (d1->start+d1->end)/2 < (d2->start+d2->end)/2;
	});

	std::vector<hts_pos_t> midpoints;
	std::vector<char*> mid_regions;
	for (deletion_t* deletion : deletions) {
		hts_pos_t midpoint = (deletion->start+deletion->end)/2;
		midpoints.push_back(midpoint);

		std::stringstream ss;
		ss << contig_name << ":" << midpoint-config.max_is << "-" << midpoint;
		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		mid_regions.push_back(region);
	}

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, mid_regions.data(), mid_regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < deletions.size() && midpoints[curr_pos] < read->core.pos) curr_pos++;

		if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
		if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0) continue;
//		if (read->core.qual < 20) continue;

		hts_pos_t start = read->core.pos + read->core.l_qseq/2;
		hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;

		for (int i = curr_pos; i < deletions.size() && midpoints[i] < read->core.pos+read->core.isize; i++) {
			if (start <= midpoints[i] && midpoints[i] <= end && read->core.isize <= config.max_is) deletions[i]->conc_pairs++;
		}
	}

	calculate_cluster_region_disc(contig_name, deletions, bam_file);

	for (char* region : mid_regions) {
		delete[] region;
	}
}


#endif //SURVINDEL2_STAT_TESTS_H
