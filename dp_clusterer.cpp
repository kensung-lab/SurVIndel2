#include <fstream>
#include <string>
#include <mutex>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "libs/cptl_stl.h"
#include "libs/IntervalTree.h"
#include "utils.h"
#include "sam_utils.h"
#include "stat_tests.h"
#include "remapping.h"
#include "libs/kdtree.h"

std::string workdir;
std::mutex mtx;

config_t config;
stats_t stats;
chr_seqs_map_t chr_seqs;
std::string reference_fname;

bcf_hdr_t* dp_vcf_header;

std::unordered_map<std::string, std::vector<deletion_t*>> deletions_by_chr;
std::unordered_map<std::string, std::vector<duplication_t*>> duplications_by_chr;
std::unordered_map<std::string, std::vector<bcf1_t*>> sr_entries_by_chr, sr_nonpass_entries_by_chr;
std::unordered_map<std::string, std::vector<bcf1_t*>> dp_entries_by_chr;
std::mutex maps_mtx;

int max_cluster_size, min_cluster_size;
std::vector<double> global_crossing_isize_dist;
std::vector<uint32_t> median_crossing_count_geqi_by_isize;

enum cluster_type_enum { LONG_PAIR, OUTWARD_PAIR };

struct cluster_t {
	hts_pos_t la_start, la_end, ra_start, ra_end;
	bool la_end_clipped, ra_start_clipped;
	int count = 1;
	uint8_t max_mapq;
	bool used = false;
	cluster_type_enum type;
	std::string rightmost_lseq, leftmost_rseq; // right-most right-facing sequence and left-most left-facing seq
	std::vector<std::string> reads;

	cluster_t() : la_start(0), la_end(0), ra_start(0), ra_end(0), max_mapq(0) {}
	cluster_t(bam1_t* read) : la_start(read->core.pos), la_end(bam_endpos(read)), ra_start(read->core.mpos), ra_end(get_mate_endpos(read)),
			max_mapq(std::max(read->core.qual, (uint8_t) get_mq(read))) {
		la_end_clipped = is_right_clipped(read, config.min_clip_len);
		ra_start_clipped = is_mate_left_clipped(read, config.min_clip_len);
		rightmost_lseq = get_sequence(read);
		leftmost_rseq = bam_get_qname(read);
		type = read->core.isize < 0 ? OUTWARD_PAIR : LONG_PAIR;
	}
};
bool operator < (const cluster_t& c1, const cluster_t& c2) {
	return std::make_tuple(c1.la_start, c1.la_end, c1.ra_start, c1.ra_end) < std::make_tuple(c2.la_start, c2.la_end, c2.ra_start, c2.ra_end);
}
bool operator == (const cluster_t& c1, const cluster_t& c2) {
	return c1.la_start == c2.la_start && c1.la_end == c2.la_end && c1.ra_start == c2.ra_start && c1.ra_end == c2.ra_end;
}
bool operator != (const cluster_t& c1, const cluster_t& c2) {
	return !(c1 == c2);
}
hts_pos_t distance(cluster_t* c1, cluster_t* c2) {
	hts_pos_t la_dist = std::max(c1->la_end, c2->la_end) - std::min(c1->la_start, c2->la_start);
	hts_pos_t ra_dist = std::max(c1->ra_end, c2->ra_end) - std::min(c1->ra_start, c2->ra_start);
	return std::max(la_dist, ra_dist);
}
cluster_t* merge(cluster_t* c1, cluster_t* c2) {
	cluster_t* merged = new cluster_t;
	merged->la_start = std::min(c1->la_start, c2->la_start);
	merged->la_end = std::max(c1->la_end, c2->la_end);
	merged->ra_start = std::min(c1->ra_start, c2->ra_start);
	merged->ra_end = std::max(c1->ra_end, c2->ra_end);
	merged->count = c1->count + c2->count;
	merged->max_mapq = std::max(c1->max_mapq, c2->max_mapq);
	merged->type = c1->type;
	merged->reads.insert(merged->reads.end(), c1->reads.begin(), c1->reads.end());
	merged->reads.insert(merged->reads.end(), c2->reads.begin(), c2->reads.end());

	// set rightmost_lseq
	if (c1->la_end > c2->la_end) merged->rightmost_lseq = c1->rightmost_lseq;
	else if (c1->la_end < c2->la_end) merged->rightmost_lseq = c2->rightmost_lseq;
	else if (c1->la_start >= c2->la_start) merged->rightmost_lseq = c1->rightmost_lseq; // both reads may be clipped at the same position, then we choose the one that starts later
	else merged->rightmost_lseq = c2->rightmost_lseq;

	// set leftmost_lseq
	if (c1->ra_start < c2->ra_start) merged->leftmost_rseq = c1->leftmost_rseq;
	else if (c1->ra_start > c2->ra_start) merged->leftmost_rseq = c2->leftmost_rseq;
	else if (c1->la_end < c2->la_end) merged->leftmost_rseq = c1->leftmost_rseq;
	else merged->leftmost_rseq = c2->leftmost_rseq;

	return merged;
}

struct cc_pair {
	cluster_t* c1, * c2;
	hts_pos_t dist;

	cc_pair(cluster_t* c1, cluster_t* c2) : c1(c1), c2(c2), dist(distance(c1, c2)) {}

	bool compatible() {
//		if (c1->la_end_clipped && c2->la_end > c1->la_end+config.min_clip_len) return false;
//		if (c2->la_end_clipped && c1->la_end > c2->la_end+config.min_clip_len) return false;
		return dist <= config.max_is;
	}
};
bool operator < (const cc_pair& cc1, const cc_pair& cc2) { return cc1.dist > cc2.dist; }

void cluster_clusters(std::vector<cluster_t*>& clusters) {
	std::sort(clusters.begin(), clusters.end(), [](cluster_t* c1, cluster_t* c2) {
		return std::make_tuple(c1->la_start, c1->la_end, c1->ra_start, c1->ra_end) < std::make_tuple(c2->la_start, c2->la_end, c2->ra_start, c2->ra_end);
	});

	std::vector<cluster_t*> dedup_clusters;
	dedup_clusters.push_back(clusters[0]);
	for (int i = 1; i < clusters.size(); i++) {
		if (*clusters[i-1] != *clusters[i]) {
			dedup_clusters.push_back(clusters[i]);
		}
	}
	clusters.swap(dedup_clusters);

	std::priority_queue<cc_pair> pq;
	for (int i = 0; i < clusters.size(); i++) {
		if (clusters[i]->used) continue;

		std::vector<cc_pair> ccps;
		for (int j = i+1; j < clusters.size(); j++) {
			if (clusters[j]->la_start-clusters[i]->la_start > config.max_is || ccps.size() > max_cluster_size) break;
			cc_pair ccp = cc_pair(clusters[i], clusters[j]);
			if (ccp.compatible()) ccps.push_back(ccp);
		}
		if (ccps.size() <= max_cluster_size) {
			for (cc_pair& ccp : ccps) pq.push(ccp);
		} else { // very large cluster - probably abnormal region - mark all involved pairs to be discarded
			clusters[i]->used = true;
			for (cc_pair& ccp : ccps) {
				ccp.c2->used = true;
			}
		}
	}

	kdtree* kd_tree_endpoints = kd_create(2);
	for (int i = 0; i < clusters.size(); i++) {
		if (clusters[i]->used) continue;
		double p[2] = {double(clusters[i]->la_end), double(clusters[i]->ra_end)};
		kd_insert(kd_tree_endpoints, p, clusters[i]);
	}

	while (!pq.empty()) {
		cc_pair ccp = pq.top();
		pq.pop();

		if (ccp.c1->used || ccp.c2->used) continue;
		ccp.c1->used = ccp.c2->used = true;

		cluster_t* merged = merge(ccp.c1, ccp.c2);
		double ps[2] = {double(merged->la_start), double(merged->ra_start)};
		kdres* res = kd_nearest_range(kd_tree_endpoints, ps, 2*config.max_is);
		while (!kd_res_end(res)) {
			cluster_t* c = (cluster_t*) kd_res_item_data(res);
			if (!c->used) {
				cc_pair ccp = cc_pair(merged, c);
				if (ccp.compatible()) pq.push(ccp);
			}
			kd_res_next(res);
		}
		kd_res_free(res);

		double pe[2] = {double(merged->la_end), double(merged->ra_end)};
		kd_insert(kd_tree_endpoints, pe, merged);
		clusters.push_back(merged);
	}
	kd_free(kd_tree_endpoints);

	std::sort(clusters.begin(), clusters.end());
}

void cluster_dps(int id, int contig_id, std::string contig_name) {
	std::string dp_fname = workdir + "/workspace/" + std::to_string(contig_id) + "-DP.bam";
	open_samFile_t* dp_bam_file = open_samFile(dp_fname, true);

	std::vector<cluster_t*> lp_clusters, ow_clusters;

	hts_itr_t* iter = sam_itr_querys(dp_bam_file->idx, dp_bam_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(dp_bam_file->file, iter, read) >= 0) {
		cluster_t* cluster = new cluster_t(read);
		cluster->reads.push_back(bam_get_qname(read));
		if (cluster->type == LONG_PAIR) {
			lp_clusters.push_back(cluster);
		} else {
//			ow_clusters.push_back(cluster);
		}
	}
	close_samFile(dp_bam_file);

	if (!lp_clusters.empty()) cluster_clusters(lp_clusters);
	if (!ow_clusters.empty()) cluster_clusters(ow_clusters);

	std::cout << "Clustered " << contig_name << std::endl;

	auto set_indel_info = [&contig_name](indel_t* indel, cluster_t* c) {
		indel->disc_pairs = c->count;
		indel->extra_info += "LC=" + std::to_string(c->la_start) + "-" + std::to_string(c->la_end) + ",";
		indel->extra_info += "RC=" + std::to_string(c->ra_start) + "-" + std::to_string(c->ra_end) + ",";

		indel->rightmost_rightfacing_seq = c->rightmost_lseq;
		indel->leftmost_leftfacing_seq = c->leftmost_rseq;
	};

	std::vector<deletion_t*> deletions;
	for (cluster_t* c : lp_clusters) {
		if (c->used || c->count < min_cluster_size || c->max_mapq != 60) continue;

		if (c->la_end < c->ra_start) {
			deletion_t* del = new deletion_t(c->la_end, c->ra_start, c->la_start, c->ra_end, NULL, NULL, 0, 0, "DP", "");
			set_indel_info(del, c);
			deletions.push_back(del);
		}
	}

	std::vector<duplication_t*> duplications;
	for (cluster_t* c : ow_clusters) {
		if (c->used || c->count < min_cluster_size || c->max_mapq != 60) continue;

		if (c->la_end > c->ra_start) {
			duplication_t* dup = new duplication_t(c->ra_start, c->la_end, c->la_start, c->ra_end, NULL, NULL, "DP");
			set_indel_info(dup, c);
			duplications.push_back(dup);
		}
	}

	maps_mtx.lock();
	deletions_by_chr[contig_name] = deletions;
	duplications_by_chr[contig_name] = duplications;
	maps_mtx.unlock();

	for (cluster_t* c : lp_clusters) {
		delete c;
	}
	for (cluster_t* c : ow_clusters) {
		delete c;
	}
}

void merge_sr_dp(int id, int contig_id, std::string contig_name, bcf_hdr_t* sr_hdr) {
	maps_mtx.lock();
	std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
	std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
	std::vector<bcf1_t*>& sr_entries = sr_entries_by_chr[contig_name];
	maps_mtx.unlock();

	if (deletions.empty() || sr_entries.empty()) return;

	std::vector<Interval<bcf1_t*> > del_intervals, dup_intervals;
	for (bcf1_t* sr_entry : sr_entries) {
		std::string sv_type = get_sv_type(sr_hdr, sr_entry);
		if (sv_type == "DEL") {
			del_intervals.push_back(Interval<bcf1_t*>(sr_entry->pos, get_sv_end(sr_hdr, sr_entry), sr_entry));
		} else {
			dup_intervals.push_back(Interval<bcf1_t*>(sr_entry->pos, get_sv_end(sr_hdr, sr_entry), sr_entry));
		}
	}
	IntervalTree<bcf1_t*> del_tree = IntervalTree<bcf1_t*>(del_intervals);
	IntervalTree<bcf1_t*> dup_tree = IntervalTree<bcf1_t*>(dup_intervals);

	for (int i = 0; i < deletions.size(); i++) {
		deletion_t* del = deletions[i];
		std::vector<Interval<bcf1_t*> > iv = del_tree.findContained(del->rc_anchor_start, del->lc_anchor_end);
		if (iv.size() != 1) continue;

		bcf1_t* corr_sr_del = iv[0].value;
		if (corr_sr_del->pos-del->rc_anchor_start <= config.max_is && del->lc_anchor_end-get_sv_end(sr_hdr, corr_sr_del) <= config.max_is) {
			if (bcf_has_filter(sr_hdr, corr_sr_del, (char*) "PASS")) {
				bcf_update_info_int32(sr_hdr, corr_sr_del, "DISC_PAIRS", &del->disc_pairs, 1);
				deletions[i] = NULL;
			} else { // if the SR deletion failed the filters, just borrow the coordinates (and the split reads)
				int* clipped_reads = NULL, n = 0;
				bcf_get_info_int32(sr_hdr, corr_sr_del, "CLIPPED_READS", &clipped_reads, &n);
				if (clipped_reads[0] == 0 || clipped_reads[1] == 0) continue; // only use 2SR
				deletions[i]->start = corr_sr_del->pos;
				deletions[i]->end = get_sv_end(sr_hdr, corr_sr_del);
				// we are only interested in the number of split reads
				deletions[i]->rc_consensus = new consensus_t(false, 0, 0, 0, 0, 0, 0, "", clipped_reads[0], 0, 0);
				deletions[i]->lc_consensus = new consensus_t(false, 0, 0, 0, 0, 0, 0, "", clipped_reads[1], 0, 0);
			}
		}
	}
	deletions.erase(std::remove(deletions.begin(), deletions.end(), (deletion_t*) NULL), deletions.end());

	for (int i = 0; i < duplications.size(); i++) {
		duplication_t* dup = duplications[i];
		std::vector<Interval<bcf1_t*> > iv = dup_tree.findContained(dup->lc_anchor_end-config.max_is,  dup->rc_anchor_start+config.max_is);
		if (iv.size() != 1) continue;

		bcf1_t* corr_sr_dup = iv[0].value;
		if (corr_sr_dup->pos < dup->lc_anchor_end && get_sv_end(sr_hdr, corr_sr_dup) > dup->rc_anchor_start) {
			if (bcf_has_filter(sr_hdr, corr_sr_dup, (char*) "PASS")) {
//				bcf_update_info_int32(sr_hdr, corr_sr_dup, "DISC_PAIRS", &dup->disc_pairs, 1);
//				duplications[i] = NULL;
			}
		} else { // if the SR deletion failed the filters, just borrow the coordinates (and the split reads)
			int* clipped_reads = NULL, n = 0;
			bcf_get_info_int32(sr_hdr, corr_sr_dup, "CLIPPED_READS", &clipped_reads, &n);
			if (clipped_reads[0] == 0 || clipped_reads[1] == 0) continue; // only use 2SR
		}
	}
	duplications.erase(std::remove(duplications.begin(), duplications.end(), (duplication_t*) NULL), duplications.end());

	std::unordered_map<std::string, deletion_t*> mateseqs_to_retrieve;
	for (deletion_t* del : deletions) {
		mateseqs_to_retrieve[del->leftmost_leftfacing_seq] = del;
	}

	std::ifstream mateseqs_fin(workdir + "/workspace/" + std::to_string(contig_id) + ".mateseqs");
	std::string qname, seq;
	while (mateseqs_fin >> qname >> seq) {
		if (!mateseqs_to_retrieve.count(qname)) continue;
		mateseqs_to_retrieve[qname]->leftmost_leftfacing_seq = seq;
		mateseqs_to_retrieve.erase(qname);
	}

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, true);
	StripedSmithWaterman::Filter filter;
	for (int i = 0; i < deletions.size(); i++) {
		deletion_t* del = deletions[i];

		suffix_prefix_aln_t spa = aln_suffix_prefix(del->rightmost_rightfacing_seq, del->leftmost_leftfacing_seq, 1, -4, config.min_clip_len);
		del->overlap = spa.overlap;
		del->mismatches = spa.mismatches;
		if (del->overlap && del->mismatches == 0) {
			std::string full_junction_seq = del->rightmost_rightfacing_seq + del->leftmost_leftfacing_seq.substr(del->overlap);

			hts_pos_t ref_lh_start = del->rc_anchor_start, ref_lh_end = del->start + full_junction_seq.length();

			hts_pos_t ref_rh_start = del->end - full_junction_seq.length(), ref_rh_end = del->lc_anchor_end;

			if (has_Ns(chr_seqs.get_seq(contig_name), ref_lh_start, ref_lh_end-ref_lh_start)) continue;
			if (has_Ns(chr_seqs.get_seq(contig_name), ref_rh_start, ref_rh_end-ref_rh_start)) continue;

			indel_t* indel = remap_consensus(full_junction_seq, chr_seqs.get_seq(contig_name), ref_lh_start, ref_lh_end-ref_lh_start, ref_rh_start,
					ref_rh_end-ref_rh_start, aligner, NULL, NULL, "DP");
			if (indel == NULL || indel->indel_type() != "DEL" || indel->len() == 0) {
				deletions[i] = NULL;
				continue;
			}

			deletions[i] = (deletion_t*) indel;
			indel->disc_pairs = del->disc_pairs;
			indel->overlap = del->overlap;
			indel->mismatches = del->mismatches;
			deletions[i]->original_range = std::to_string(del->start) + "-" + std::to_string(del->end);
			deletions[i]->remapped = true;
		}
	}
	deletions.erase(std::remove(deletions.begin(), deletions.end(), (deletion_t*) NULL), deletions.end());
}

void add_filtering_info(int id, std::string contig_name, std::string bam_fname) {
	maps_mtx.lock();
	std::vector<deletion_t*>& deletions = deletions_by_chr[contig_name];
	std::vector<duplication_t*>& duplications = duplications_by_chr[contig_name];
	maps_mtx.unlock();
	if (deletions.empty()) return;

	std::vector<deletion_t*> shorter_deletions, longer_deletions;
	for (deletion_t* deletion : deletions) {
		if (deletion->len() < config.max_is) shorter_deletions.push_back(deletion);
		else longer_deletions.push_back(deletion);
	}

	open_samFile_t* bam_file = open_samFile(bam_fname, false);
	hts_set_fai_filename(bam_file->file, reference_fname.c_str());
	depth_filter_del(contig_name, deletions, bam_file, config.min_size_for_depth_filtering, config);
	depth_filter_dup(contig_name, duplications, bam_file, config.min_size_for_depth_filtering, config);
	if (!shorter_deletions.empty()) calculate_confidence_interval_size(contig_name, global_crossing_isize_dist, median_crossing_count_geqi_by_isize, shorter_deletions, bam_file, config, stats);
	if (!longer_deletions.empty()) calculate_ptn_ratio(contig_name, longer_deletions, bam_file, config);
	calculate_cluster_region_disc(contig_name, duplications, bam_file, config);
	close_samFile(bam_file);
	std::cout << "Stats calculated for " << contig_name << std::endl;

	mtx.lock();
	bcf1_t* bcf_entry = bcf_init();
	for (deletion_t* del : deletions) {
		std::vector<std::string> filters;

		if (del->len() < config.min_sv_size) {
			continue;
			filters.push_back("SMALL");
		}
		if (del->ks_pval > 0.01) {
			filters.push_back("KS_FILTER");
		} else
		if (del->len() >= config.max_is && double(del->disc_pairs)/(del->disc_pairs+del->conc_pairs) < 0.25) {
			filters.push_back("LOW_PTN_RATIO");
		}

		if (del->med_left_flanking_cov*0.74<=del->med_indel_left_cov || del->med_right_flanking_cov*0.74<=del->med_indel_right_cov) {
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

		if (del->len() > 10000 && (del->l_cluster_region_disc_pairs >= del->disc_pairs || del->r_cluster_region_disc_pairs >= del->disc_pairs)) {
			filters.push_back("AMBIGUOUS_REGION");
		}

		if (filters.empty()) {
			filters.push_back("PASS");
		}

		del2bcf(dp_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, del, filters);
		float ks_pval = del->ks_pval;
		bcf_update_info_float(dp_vcf_header, bcf_entry, "KS_PVAL", &ks_pval, 1);
		bcf_update_info_int32(dp_vcf_header, bcf_entry, "OVERLAP", &del->overlap, 1);
		bcf_update_info_int32(dp_vcf_header, bcf_entry, "MISMATCHES", &del->mismatches, 1);
//		if (!del->genotype.empty()) bcf_update_info_string(dp_vcf_header, bcf_entry, "GENOTYPE", del->genotype.c_str());
		if (!del->original_range.empty()) bcf_update_info_string(dp_vcf_header, bcf_entry, "ORIGINAL_RANGE", del->original_range.c_str());
		dp_entries_by_chr[contig_name].push_back(bcf_dup(bcf_entry));
	}

	for (duplication_t* dup : duplications) {
		std::vector<std::string> filters;

		if (dup->med_left_flanking_cov*1.25>dup->med_indel_left_cov || dup->med_right_flanking_cov*1.25>dup->med_indel_right_cov) {
			filters.push_back("DEPTH_FILTER");
		}
		if (dup->med_left_flanking_cov > stats.max_depth || dup->med_right_flanking_cov > stats.max_depth ||
			dup->med_left_flanking_cov < stats.min_depth || dup->med_right_flanking_cov < stats.min_depth) {
			filters.push_back("ANOMALOUS_FLANKING_DEPTH");
		}

		if (dup->lc_cluster_region_disc_pairs >= dup->disc_pairs || dup->rc_cluster_region_disc_pairs >= dup->disc_pairs) {
			filters.push_back("AMBIGUOUS_REGION");
		}

		if (filters.empty()) {
			filters.push_back("PASS");
		}

		dup2bcf(dp_vcf_header, bcf_entry, chr_seqs.get_seq(contig_name), contig_name, dup, filters);
		dp_entries_by_chr[contig_name].push_back(bcf_dup(bcf_entry));
	}

	bcf_destroy(bcf_entry);
	mtx.unlock();
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];

    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    reference_fname = argv[3];
    std::string sample_name = argv[4];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");

	stats.parse(workdir + "/stats.txt");
	min_cluster_size = std::max(3, int(stats.avg_depth+5)/10);

	max_cluster_size = (stats.max_depth * config.max_is)/config.read_len;

	chr_seqs.read_fasta_into_map(reference_fname);

	std::ifstream crossing_isizes_dist_fin(workdir + "/crossing_isizes.txt");
	int isize, count;
	while (crossing_isizes_dist_fin >> isize >> count) {
		for (int i = 0; i < count; i++) global_crossing_isize_dist.push_back(isize);
	}
	std::random_shuffle(global_crossing_isize_dist.begin(), global_crossing_isize_dist.end());
	global_crossing_isize_dist.resize(100000);
	crossing_isizes_dist_fin.close();

	std::ifstream crossing_isizes_count_geq_i_fin(workdir + "/crossing_isizes_count_geq_i.txt");
	int median;
	while (crossing_isizes_count_geq_i_fin >> isize >> median) {
		median_crossing_count_geqi_by_isize.push_back(median);
	}
	crossing_isizes_count_geq_i_fin.close();

	dp_vcf_header = generate_vcf_header(chr_seqs, sample_name, config, full_cmd_str);
	std::string dp_vcf_fname = workdir + "/dp.vcf.gz";
	htsFile* dp_vcf_file = bcf_open(dp_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(dp_vcf_file, dp_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + dp_vcf_fname + ".");
	}

	ctpl::thread_pool thread_pool1(config.threads);
	std::vector<std::future<void> > futures;
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool1.push(cluster_dps, contig_id, contig_name);
		futures.push_back(std::move(future));
	}
	thread_pool1.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		try {
			futures[i].get();
		} catch (char const* s) {
			std::cout << s << std::endl;
		}
	}
	futures.clear();

	std::string sr_vcf_fname = workdir + "/sr.norm.dedup.vcf.gz";
	htsFile* sr_vcf_file = bcf_open(sr_vcf_fname.c_str(), "r");
	bcf_hdr_t* sr_vcf_hdr = bcf_hdr_read(sr_vcf_file);
	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(sr_vcf_file, sr_vcf_hdr, bcf_entry) == 0) {
		int* clipped_reads = NULL, n = 0;
		bcf_get_info_int32(sr_vcf_hdr, bcf_entry, "CLIPPED_READS", &clipped_reads, &n);
		if (bcf_has_filter(sr_vcf_hdr, bcf_entry, (char*) "PASS") || (clipped_reads[0] == 0 && clipped_reads[1] == 0)) { // only use 2SR or PASS
			sr_entries_by_chr[bcf_seqname_safe(sr_vcf_hdr, bcf_entry)].push_back(bcf_dup(bcf_entry));
		} else {
			sr_nonpass_entries_by_chr[bcf_seqname_safe(sr_vcf_hdr, bcf_entry)].push_back(bcf_dup(bcf_entry));
		}
	}
	bcf_close(sr_vcf_file);

	// merge SR and DISC
	ctpl::thread_pool thread_pool2(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool2.push(merge_sr_dp, contig_id, contig_name, sr_vcf_hdr);
		futures.push_back(std::move(future));
	}
	thread_pool2.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		try {
			futures[i].get();
		} catch (char const* s) {
			std::cout << s << std::endl;
		}
	}
	futures.clear();

	ctpl::thread_pool thread_pool3(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool3.push(add_filtering_info, contig_name, bam_fname);
		futures.push_back(std::move(future));
	}
	thread_pool3.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		try {
			futures[i].get();
		} catch (char const* s) {
			std::cout << s << std::endl;
		}
	}

	int del_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
    	auto& bcf_entries_contig = dp_entries_by_chr[contig_name];
    	std::sort(bcf_entries_contig.begin(), bcf_entries_contig.end(), [](const bcf1_t* b1, const bcf1_t* b2) {return b1->pos < b2->pos;});
		for (bcf1_t* bcf_entry : bcf_entries_contig) {
			std::string id = "DEL_DP_" + std::to_string(del_id++);
			bcf_update_id(dp_vcf_header, bcf_entry, id.c_str());
			if (bcf_write(dp_vcf_file, dp_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + dp_vcf_fname + ".");
			}
		}
    }

	bcf_close(dp_vcf_file);

	tbx_index_build(dp_vcf_fname.c_str(), 0, &tbx_conf_vcf);


	std::string merged_vcf_fname = workdir + "/out.vcf.gz";
	htsFile* merged_vcf_file = bcf_open(merged_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(merged_vcf_file, dp_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + merged_vcf_fname + ".");
	}

	std::string merged_pass_vcf_fname = workdir + "/out.pass.vcf.gz";
	htsFile* merged_pass_vcf_file = bcf_open(merged_pass_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(merged_pass_vcf_file, dp_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + merged_pass_vcf_fname + ".");
	}

	for (std::string& contig_name : chr_seqs.ordered_contigs) {
		std::vector<bcf1_t*> all_entries;
		all_entries.insert(all_entries.end(), dp_entries_by_chr[contig_name].begin(), dp_entries_by_chr[contig_name].end());
		all_entries.insert(all_entries.end(), sr_entries_by_chr[contig_name].begin(), sr_entries_by_chr[contig_name].end());
		all_entries.insert(all_entries.end(), sr_nonpass_entries_by_chr[contig_name].begin(), sr_nonpass_entries_by_chr[contig_name].end());
		std::sort(all_entries.begin(), all_entries.end(), [](bcf1_t* b1, bcf1_t* b2) {
			return b1->pos < b2->pos;
		});
		for (bcf1_t* bcf_entry : all_entries) {
			if (bcf_write(merged_vcf_file, dp_vcf_header, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + merged_vcf_fname + ".");
			}
			if (bcf_has_filter(dp_vcf_header, bcf_entry, (char*) "PASS")) {
				if (bcf_write(merged_pass_vcf_file, dp_vcf_header, bcf_entry) != 0) {
					throw std::runtime_error("Failed to write to " + merged_vcf_fname + ".");
				}
			}
		}
	}

	bcf_close(merged_vcf_file);
	bcf_close(merged_pass_vcf_file);

	tbx_index_build(merged_vcf_fname.c_str(), 0, &tbx_conf_vcf);
	tbx_index_build(merged_pass_vcf_fname.c_str(), 0, &tbx_conf_vcf);

}
