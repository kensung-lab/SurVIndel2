#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "sam_utils.h"
#include "utils.h"

void merge_hsr_sr_with_1sr_rc(bcf_hdr_t* hdr, bcf1_t* hsr_b, bcf1_t* sr_b, bool is_sr_rc) {
	int* hsr_cr = NULL;
	int size = 0;
	bcf_get_info_int32(hdr, hsr_b, "CLIPPED_READS", &hsr_cr, &size);

	int* sr_cr = NULL;
	size = 0;
	bcf_get_info_int32(hdr, sr_b, "CLIPPED_READS", &sr_cr, &size);

	int* hsr_mmapq = NULL;
	size = 0;
	bcf_get_info_int32(hdr, hsr_b, "MAX_MAPQ", &hsr_mmapq, &size);

	int* sr_mmapq = NULL;
	size = 0;
	bcf_get_info_int32(hdr, sr_b, "MAX_MAPQ", &sr_mmapq, &size);

	std::string sv_type = get_sv_type(hdr, sr_b);
	if (sv_type == "DEL") { // RC = 1st bkp, LC = 2nd bkp
		int idx = is_sr_rc ? 0 : 1;
		hsr_cr[idx] = sr_cr[idx];
		hsr_mmapq[idx] = sr_mmapq[idx];
	} else if (sv_type == "DUP") { // RC = 2nd bkp, LC = 1st bkp
		int idx = is_sr_rc ? 1 : 0;
		hsr_cr[idx] = sr_cr[idx];
		hsr_mmapq[idx] = sr_mmapq[idx];
	}

	bcf_update_info_int32(hdr, hsr_b, "CLIPPED_READS", hsr_cr, 2);
	bcf_update_info_int32(hdr, hsr_b, "MAX_MAPQ", hsr_mmapq, 2);
	bcf_update_info_string(hdr, hsr_b, "SOURCE", "2SR");
}

void merge_sr_rc_with_sr_lc(bcf_hdr_t* hdr, bcf1_t* rc_b, bcf1_t* lc_b) {
	merge_hsr_sr_with_1sr_rc(hdr, rc_b, lc_b, false);

	bcf_update_info_int32(hdr, rc_b, "WEAK_SUPPORTING_1SR_READS", NULL, 0);
	bcf_update_info_int32(hdr, rc_b, "STRONG_SUPPORTING_1SR_READS", NULL, 0);
	bcf_update_info_int32(hdr, rc_b, "PERFECT_SUPPORTING_1SR_READS", NULL, 0);
}

void merge_sr_lc_with_sr_rc(bcf_hdr_t* hdr, bcf1_t* lc_b, bcf1_t* rc_b) {
	merge_hsr_sr_with_1sr_rc(hdr, lc_b, rc_b, true);

	bcf_update_info_int32(hdr, lc_b, "WEAK_SUPPORTING_1SR_READS", NULL, 0);
	bcf_update_info_int32(hdr, lc_b, "STRONG_SUPPORTING_1SR_READS", NULL, 0);
	bcf_update_info_int32(hdr, lc_b, "PERFECT_SUPPORTING_1SR_READS", NULL, 0);
}

int main(int argc, char* argv[]) {

	std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];

	std::unordered_map<std::string, int> source_priorities;
	source_priorities["2SR"] = 1;
	source_priorities["HSR-SR"] = 2;
	source_priorities["SR-HSR"] = 2;
	source_priorities["2HSR"] = 3;
	source_priorities["1SR_LC"] = 4;
	source_priorities["1SR_RC"] = 4;
	source_priorities["1HSR_LC"] = 5;
	source_priorities["1HSR_RC"] = 5;
	source_priorities["DP"] = 6;

	std::unordered_map<std::string, bcf1_t*> bcf_entries;
	std::vector<std::string> ids_sorted;

	htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	bcf_hdr_t* in_vcf_hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* b = bcf_init();
	while (bcf_read(in_vcf_file, in_vcf_hdr, b) == 0) {
		std::string unique_id = std::string(bcf_seqname_safe(in_vcf_hdr, b)) + "_" + get_sv_type(in_vcf_hdr, b) + "_" +
				std::to_string(b->pos) + "_" + std::to_string(get_sv_end(in_vcf_hdr, b));
		if (bcf_entries.count(unique_id)) {
			// if this is not PASS and there is an equivalent indel already stored, ignore
			if (!bcf_has_filter(in_vcf_hdr, b, (char*) "PASS")) continue;
			else if (!bcf_has_filter(in_vcf_hdr, bcf_entries[unique_id], (char*) "PASS")) {
				// if stored indel is not PASS, but there is an equivalent PASS indel, store that instead
				bcf_entries[unique_id] = bcf_dup(b);
				continue;
			}

			std::string src_prev = get_sv_info_str(in_vcf_hdr, bcf_entries[unique_id], "SOURCE");
			std::string src_curr = get_sv_info_str(in_vcf_hdr, b, "SOURCE");
			int src_prev_pri = source_priorities[src_prev];
			int src_curr_pri = source_priorities[src_curr];

			if (src_prev_pri > src_curr_pri) {
				std::swap(bcf_entries[unique_id], b);
				std::swap(src_prev, src_curr);
				std::swap(src_prev_pri, src_curr_pri);
			}

			if (src_prev == "2SR") continue; // 2SR has highest priority. If prev indel is 2SR, keep it
			else if (src_prev == "HSR-SR" && src_curr == "1SR_RC") { // HSR-SR + 1SR_RC -> 2SR
				merge_hsr_sr_with_1sr_rc(in_vcf_hdr, bcf_entries[unique_id], b, true);
			} else if (src_prev == "SR-HSR" && src_curr == "1SR_LC") { // SR-HSR + 1SR_LC -> 2SR
				merge_hsr_sr_with_1sr_rc(in_vcf_hdr, bcf_entries[unique_id], b, false);
			} else if (src_prev == "1SR_LC" && src_curr == "1SR_RC") { // 1SR_RC + 1SR_LC -> 2SR
				merge_sr_lc_with_sr_rc(in_vcf_hdr, bcf_entries[unique_id], b);
			} else if (src_prev == "1SR_RC" && src_curr == "1SR_LC") { // 1SR_RC + 1SR_LC -> 2SR
				merge_sr_rc_with_sr_lc(in_vcf_hdr, bcf_entries[unique_id], b);
			} else if (src_prev == "1SR_RC" && src_curr == "1HSR_LC") { // 1SR_RC + 1HSR_LC -> SR-HSR
				merge_sr_rc_with_sr_lc(in_vcf_hdr, bcf_entries[unique_id], b);
				bcf_update_info_string(in_vcf_hdr, bcf_entries[unique_id], "SOURCE", "SR-HSR");
			} else if (src_prev == "1SR_LC" && src_curr == "1HSR_RC") { // 1HSR_RC + 1SR_LC -> HSR-SR
				merge_sr_lc_with_sr_rc(in_vcf_hdr, bcf_entries[unique_id], b);
				bcf_update_info_string(in_vcf_hdr, bcf_entries[unique_id], "SOURCE", "HSR-SR");
			} else if (src_prev == "1HSR_RC" && src_curr == "1HSR_LC") { // 1HSR_RC + 1HSR_LC -> 2HSR
				merge_sr_rc_with_sr_lc(in_vcf_hdr, bcf_entries[unique_id], b);
				bcf_update_info_string(in_vcf_hdr, bcf_entries[unique_id], "SOURCE", "2HSR");
			} else if (src_prev == "1HSR_LC" && src_curr == "1HSR_RC") { // 1HSR_RC + 1HSR_LC -> 2HSR
				merge_sr_lc_with_sr_rc(in_vcf_hdr, bcf_entries[unique_id], b);
				bcf_update_info_string(in_vcf_hdr, bcf_entries[unique_id], "SOURCE", "2HSR");
			}

			// TODO: other possible cases: HSR-SR + SR-HSR -> 2SR, 2HSR + 1SR -> HSR-SR (or SR-HSR)
		} else {
			ids_sorted.push_back(unique_id);
			bcf_entries[unique_id] = bcf_dup(b);
		}
	}

	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, in_vcf_hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	for (std::string& id : ids_sorted) {
		if (bcf_write(out_vcf_file, in_vcf_hdr, bcf_entries[id]) != 0) {
			throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
		}
	}

	bcf_close(out_vcf_file);
}
