#ifndef SURVINDEL2_UTILS_H
#define SURVINDEL2_UTILS_H

#include <atomic>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "htslib/kseq.h"
#include <unistd.h>
#include <numeric>
#include <chrono>
#include <ctime>
KSEQ_INIT(int, read)

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "libs/IntervalTree.h"
#include "libs/ssw_cpp.h"
#include "libs/ssw.h"

struct consensus_t {
    bool left_clipped;
    int contig_id;
    hts_pos_t breakpoint, start, end; // we follow the vcf conventions, i.e. this is the base "before" the breakpoint
    int clip_len, lowq_clip_portion;
    std::string consensus;
    size_t supp_clipped_reads;
    uint8_t max_mapq;
    hts_pos_t remap_boundary;
    bool is_hsr = false;

    static const int LOWER_BOUNDARY_NON_CALCULATED = 0, UPPER_BOUNDARY_NON_CALCULATED = INT32_MAX;
    static const int UNKNOWN_CLIP_LEN = INT16_MAX;

    consensus_t(bool left_clipped, int contig_id, hts_pos_t start, hts_pos_t end, hts_pos_t breakpoint, int clip_len, int lowq_clip_portion,
                const std::string& consensus, size_t supp_clipped_reads, uint8_t max_mapq, hts_pos_t remap_boundary)
            : left_clipped(left_clipped), contig_id(contig_id), start(start), end(end), breakpoint(breakpoint), clip_len(clip_len),
              lowq_clip_portion(lowq_clip_portion), consensus(consensus), supp_clipped_reads(supp_clipped_reads),
			  max_mapq(max_mapq), remap_boundary(remap_boundary) {}

    char dir() { return left_clipped ? 'L' : 'R'; }

    std::string name() {
        std::stringstream ss;
        ss << contig_id << "_" << breakpoint << "_" << dir() << "_" << clip_len << "_" << lowq_clip_portion << "_" << supp_clipped_reads;
        return ss.str();
    }

    int anchor_len() { return consensus.length() - clip_len; }
};

std::atomic<int> _indel_id;
struct sr_remap_info_t {
	hts_pos_t remapped_bp; // breakpoint that was found by remapping the clip
	hts_pos_t remapped_other_end;
	int weak_supporting_reads = 0, strong_supporting_reads = 0, perfect_supporting_reads = 0, tried_reads = 0;

	sr_remap_info_t(hts_pos_t remapped_bp, hts_pos_t remapped_other_end) : remapped_bp(remapped_bp), remapped_other_end(remapped_other_end) {}
};
struct indel_t {
    std::string id;
    hts_pos_t start, end;
    hts_pos_t rc_anchor_start, lc_anchor_end; // start of the rc anchor and lc anchor end
    int disc_pairs = 0;
    consensus_t* lc_consensus,* rc_consensus;
    double mm_rate = 0.0;
    std::string source;
    int left_flanking_cov = 0, indel_left_cov = 0, indel_right_cov = 0, right_flanking_cov = 0;
    int med_left_flanking_cov = 0, med_indel_left_cov = 0, med_indel_right_cov = 0, med_right_flanking_cov = 0;
    int med_left_cluster_cov = 0, med_right_cluster_cov = 0;
    int full_junction_score = 0, split_junction_score = 0;
    std::string extra_info;
    int overlap = 0, mismatches = 0;
    bool remapped = false;
    std::string rightmost_rightfacing_seq, leftmost_leftfacing_seq;

    sr_remap_info_t* sr_remap_info = NULL;

    indel_t(hts_pos_t start, hts_pos_t end, hts_pos_t rc_anchor_start, hts_pos_t lc_anchor_end, std::string source,
            consensus_t* lc_consensus, consensus_t* rc_consensus)
    : start(start), end(end), rc_anchor_start(rc_anchor_start), lc_anchor_end(lc_anchor_end),
	  source(source), lc_consensus(lc_consensus), rc_consensus(rc_consensus) { }

    virtual hts_pos_t len() { return end-start; }
    virtual std::string indel_type() { return ""; };
    virtual std::string to_string(std::string contig_name) {
    	return contig_name + ":" + std::to_string(start) + "-" + std::to_string(end);
    }

    bool is_single_consensus() { return lc_consensus == NULL || rc_consensus == NULL; }
    bool imprecise() { return lc_consensus == NULL && rc_consensus == NULL && remapped == false; }
};

struct deletion_t : indel_t {
	static const int SIZE_NOT_COMPUTED = INT32_MAX;
	static const int REMAP_LB_NOT_COMPUTED = 0, REMAP_UB_NOT_COMPUTED = INT32_MAX;

    int lh_score, rh_score;
    int max_conf_size = SIZE_NOT_COMPUTED, estimated_size = SIZE_NOT_COMPUTED, conc_pairs = 0;
    double ks_pval = -1.0;
    hts_pos_t remap_boundary_lower = REMAP_LB_NOT_COMPUTED, remap_boundary_upper = REMAP_UB_NOT_COMPUTED;
    std::string ins_seq;
    int l_cluster_region_disc_pairs = 0, r_cluster_region_disc_pairs = 0;

    std::string original_range;
    std::string genotype;

    deletion_t(hts_pos_t start, hts_pos_t end, hts_pos_t rc_anchor_start, hts_pos_t lc_anchor_end, consensus_t* lc_consensus,
               consensus_t* rc_consensus, int lh_score, int rh_score, std::string source, std::string ins_seq) :
            indel_t(start, end, rc_anchor_start, lc_anchor_end, source, lc_consensus, rc_consensus), lh_score(lh_score), rh_score(rh_score),
			ins_seq(ins_seq) {};

    hts_pos_t len() override { return end-start-ins_seq.length(); }

    std::string indel_type() override { return "DEL"; }
};

struct duplication_t : indel_t {
    static const int OW_NOT_COMPUTED = INT32_MAX;

    hts_pos_t original_start, original_end;
    int lc_cluster_region_disc_pairs = 0, rc_cluster_region_disc_pairs = 0;

    duplication_t(hts_pos_t start, hts_pos_t end, hts_pos_t rc_anchor_start, hts_pos_t lc_anchor_end, consensus_t* lc_consensus,
                  consensus_t* rc_consensus, std::string source) :
                  indel_t(std::max(start, hts_pos_t(0)), end, rc_anchor_start, lc_anchor_end, source, lc_consensus, rc_consensus) {}
    // Normally we store/report the base BEFORE the event, as per VCF specifications.
    // We handle here the exception where the duplication starts at the first base of the contig (base before 0 is -1, but it is not valid,
    // so we still store 0)

    std::string indel_type() override { return "DUP"; }
};


struct config_t {

    int threads, seed;
    int max_is; // find a way to move this to stats_t
    int min_sv_size;
    int match_score;
    int read_len; // this is not exactly "config", but it is more convenient to place it here
    int min_clip_len;
    double max_seq_error;
    int max_clipped_pos_dist;
    int min_size_for_depth_filtering;
    std::string sampling_regions, version;

    int clip_penalty = 7;
    int min_score_diff = 15;
    int high_confidence_mapq = 60;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        threads = std::stoi(config_params["threads"]);
        seed = std::stoi(config_params["seed"]);
        max_is = std::stoi(config_params["max_is"]);
        min_sv_size = std::stoi(config_params["min_sv_size"]);
        match_score = std::stoi(config_params["match_score"]);
        min_clip_len = std::stoi(config_params["min_clip_len"]);
        read_len = std::stod(config_params["read_len"]);
        max_seq_error = std::stod(config_params["max_seq_error"]);
        max_clipped_pos_dist = std::stoi(config_params["max_clipped_pos_dist"]);
        min_size_for_depth_filtering = std::stoi(config_params["min_size_for_depth_filtering"]);
        sampling_regions = config_params["sampling_regions"];
        version = config_params["version"];
    }

    std::string clip_penalty_padding() { return std::string(this->clip_penalty, 'N'); }
};

struct stats_t {
    double avg_depth, lt_depth_stddev;
    int min_depth, max_depth;
    int pop_avg_crossing_is = 0;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        avg_depth = std::stod(config_params["avg_depth"]);
        lt_depth_stddev = std::stod(config_params["lt_depth_stddev"]);
        pop_avg_crossing_is = std::stoi(config_params["pop_avg_crossing_is"]);
        min_depth = std::stoi(config_params["min_depth"]);
        max_depth = std::stoi(config_params["max_depth"]);
    }
};

struct tandem_rep_t {

        std::string chr;
        hts_pos_t start, end;

        tandem_rep_t(std::string& line) {
            std::stringstream ss(line);
            ss >> chr >> start >> end;
            start--; end--; // annotations files are in 1-based coordinates, here we are working with 0-based
        }
};

struct tandem_rep_trees_t {

    std::unordered_map<std::string, IntervalTree<tandem_rep_t>* > trees;

    tandem_rep_trees_t() {}
    tandem_rep_trees_t(std::unordered_map<std::string, IntervalTree<tandem_rep_t>* > trees) : trees(trees) {}

    bool overlaps_tr(std::string contig, hts_pos_t start, hts_pos_t end) {
        std::vector<Interval<tandem_rep_t> > reps = trees[contig]->findOverlapping(start, end);
        return !reps.empty();
    }
};

tandem_rep_trees_t parse_tr(std::string fname, int min_sv_size) {
    std::unordered_map<std::string, IntervalTree<tandem_rep_t>* > repeat_trees;
    std::ifstream rep_fin(fname);
    if (!rep_fin.is_open()) {
        throw "ERROR: Could not open simpleRepeats file";
    } else {
        std::unordered_map<std::string, std::vector<Interval<tandem_rep_t> > > v;
        std::string line;
        while (getline(rep_fin, line)) {
            if (line[0] != '#') {
                tandem_rep_t rep(line);
                if (rep.end-rep.start >= min_sv_size) {
                    v[rep.chr].push_back(Interval<tandem_rep_t>(rep.start, rep.end, rep));
                }

            }
        }

        for (auto& k : v) {
            repeat_trees[k.first] = new IntervalTree<tandem_rep_t>(k.second);
        }
    }
    return tandem_rep_trees_t(repeat_trees);
}


struct contig_map_t {

    std::unordered_map<std::string, size_t> name_to_id;
    std::vector<std::string> id_to_name;

    contig_map_t() {}
    contig_map_t(std::string workdir) { load(workdir); }

    void load(std::string workdir) {
        std::ifstream fin(workdir + "/contig_map");
        std::string name;
        int id = 0;
        while (fin >> name) {
            name_to_id[name] = id;
            id_to_name.push_back(name);
            id++;
        }
    }

    size_t size() {return id_to_name.size();}
    std::string get_name(size_t id) {return id_to_name[id];};
};


struct chr_seq_t {
    char* seq;
    hts_pos_t len;

    chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
    ~chr_seq_t() {delete[] seq;}
};
struct chr_seqs_map_t {
    std::unordered_map<std::string, chr_seq_t*> seqs;
    std::vector<std::string> ordered_contigs;

    void read_fasta_into_map(std::string& reference_fname) {
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l);
            ordered_contigs.push_back(seq_name);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    char* get_seq(std::string seq_name) {
        return seqs[seq_name]->seq;
    }

    hts_pos_t get_len(std::string seq_name) {
        return seqs[seq_name]->len;
    }

    void clear() {
        for (auto& e : seqs) {
            delete e.second;
            e.second = NULL;
        }
    }

    ~chr_seqs_map_t() {
        clear();
    }
};

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

void strrev(char* str) {
    if (str == NULL) return;

    int len = strlen(str);
    for (int i = 0; i < len/2; i++) {
        std::swap(str[i], str[len-1-i]);
    }
}

// excludes soft-clipped parts from prefix_len
int swaln_mismatches_in_prefix(StripedSmithWaterman::Alignment& aln, int prefix_len) {
    int mismatches = 0;
    for (uint32_t op : aln.cigar) {
        if (prefix_len <= 0) break;
        char opchr = cigar_int_to_op(op);
        int oplen = cigar_int_to_len(op);
        if (opchr == 'X') {
            mismatches += std::min(oplen, prefix_len);
        }
        if (opchr != 'S') {
            prefix_len -= oplen;
        }
    }
    return mismatches;
}
int swaln_mismatches_in_suffix(StripedSmithWaterman::Alignment& aln, int suffix_len) {
    int mismatches = 0;
    for (int i = aln.cigar.size()-1; i >= 0 && suffix_len > 0; i--) {
        uint32_t op = aln.cigar[i];
        char opchr = cigar_int_to_op(op);
        int oplen = cigar_int_to_len(op);
        if (opchr == 'X') {
            mismatches += std::min(oplen, suffix_len);
        }
        if (opchr != 'S') {
            suffix_len -= oplen;
        }
    }
    return mismatches;
}

int64_t overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

struct suffix_prefix_aln_t {
    int overlap, score, mismatches;

    suffix_prefix_aln_t(int overlap, int score, int mismatches) : overlap(overlap), score(score), mismatches(mismatches) {}

    double mismatch_rate() { return overlap > 0 ? double(mismatches)/overlap : 1; }
};

// Finds the best alignment between a suffix of s1 and a prefix of s2
// Disallows gaps
suffix_prefix_aln_t aln_suffix_prefix(std::string& s1, std::string& s2, int match_score, int mismatch_score,
                                      int min_overlap = 1, int max_overlap = INT32_MAX) {
    int best_score = 0, best_aln_mismatches = 0;
    int overlap = 0;

    for (int i = std::max(0, (int) s1.length()-max_overlap); i < s1.length()-min_overlap+1; i++) {
        if (i+s2.length() < s1.length()) continue;

        int sp_len = s1.length()-i;
        if (best_score >= sp_len*match_score) break; // current best score is unbeatable

        int mismatches = 0;
        const char* s1_suffix = s1.data()+i;
        const char* s2_prefix = s2.data();
        while (*s1_suffix) {
            if (*s1_suffix != *s2_prefix) mismatches++;
            s1_suffix++; s2_prefix++;
        }

        int score = (sp_len-mismatches)*match_score + mismatches*mismatch_score;

        if (best_score < score) {
            best_score = score;
            best_aln_mismatches = mismatches;
            overlap = sp_len;
        }
    }
    return suffix_prefix_aln_t(overlap, best_score, best_aln_mismatches);
}

template<typename T>
T mean(std::vector<T>& v) {
    return std::accumulate(v.begin(), v.end(), (T)0.0)/v.size();
}

bcf_hrec_t* generate_contig_hrec() {
	bcf_hrec_t* contig_hrec = new bcf_hrec_t;
	contig_hrec->type = BCF_HL_CTG;
	contig_hrec->key = strdup("contig");
	contig_hrec->value = NULL;
	contig_hrec->keys = contig_hrec->vals = NULL;
	contig_hrec->nkeys = 0;
	int r1 = bcf_hrec_add_key(contig_hrec, "ID", 2);
	int r2 = bcf_hrec_add_key(contig_hrec, "length", 6);
	if (r1 || r2) {
		throw std::runtime_error("Failed to create contig to VCF header.");
	}
	return contig_hrec;
}
bcf_hdr_t* generate_vcf_header(chr_seqs_map_t& contigs, std::string& sample_name, config_t config, std::string command) {
	bcf_hdr_t* header = bcf_hdr_init("w");

	// add contigs
	for (std::string contig_name : contigs.ordered_contigs) {
		bcf_hrec_t* hrec = generate_contig_hrec();
		int r1 = bcf_hrec_set_val(hrec, 0, contig_name.c_str(), contig_name.length(), false);
		std::string len_str = std::to_string(contigs.get_len(contig_name));
		int r2 = bcf_hrec_set_val(hrec, 1, len_str.c_str(), len_str.length(), false);
		if (r1 || r2) {
			throw std::runtime_error("Failed to create contig to VCF header.");
		}
		bcf_hdr_add_hrec(header, hrec);
	}

	int len;

	// add FILTER tags
	const char* size_flt_tag = "##FILTER=<ID=SIZE_FILTER,Description=\"Size of the event is outside the predicted confidence interval.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, size_flt_tag, &len));

	const char* remap_boundary_flt_tag = "##FILTER=<ID=REMAP_BOUNDARY_FILTER,Description=\"One of the breakpoints is incompatible with mate locations.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_boundary_flt_tag, &len));

	const char* depth_flt_tag = "##FILTER=<ID=DEPTH_FILTER,Description=\"Depth of the region is incompatible with type of the event. Only applicable to long events.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, depth_flt_tag, &len));

	const char* anom_flanking_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_FLANKING_DEPTH,Description=\"Depth of region(s) flanking this event is anomalous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anom_flanking_depth_flt_tag, &len));

	const char* anom_del_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_DEL_DEPTH,Description=\"Depth of the deleted region is anomalous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anom_del_depth_flt_tag, &len));

	const char* not_enough_ow_pairs_flt_tag = "##FILTER=<ID=NOT_ENOUGH_OW_PAIRS,Description=\"Not enough outward oriented pairs supporting a large duplication.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, not_enough_ow_pairs_flt_tag, &len));

	const char* weak_split_aln_flt_tag = "##FILTER=<ID=WEAK_SPLIT_ALIGNMENT,Description=\"Split alignment not significantly better than full junction alignment.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_split_aln_flt_tag, &len));

	const char* low_mapq_cons_flt_tag = "##FILTER=<ID=LOW_MAPQ_CONSENSUSES,Description=\"No high MAPQ read supports the consensus(es) used to call this SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_mapq_cons_flt_tag, &len));

	const char* weak_support_flt_tag = "##FILTER=<ID=WEAK_SUPPORT,Description=\"Remapped breakpoint has low support from local reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_support_flt_tag, &len));

	const char* low_ptn_ratio_flt_tag = "##FILTER=<ID=LOW_PTN_RATIO,Description=\"Low positive-to-negative ratio.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_ptn_ratio_flt_tag, &len));

	const char* ks_filter_flt_tag = "##FILTER=<ID=KS_FILTER,Description=\"According to KS test, the local IS distribution is not "
			"sufficiently different from the global IS distribution.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ks_filter_flt_tag, &len));

	const char* ambiguous_flt_tag = "##FILTER=<ID=AMBIGUOUS_REGION,Description=\"Region containing the deletion is ambiguous.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ambiguous_flt_tag, &len));

	// add INFO tags
	const char* svtype_tag = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the indel (DEL or DUP).\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svtype_tag, &len));

	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, end_tag, &len));

	const char* svlen_tag = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svlen_tag, &len));

	const char* max_size_tag = "##INFO=<ID=MAX_SIZE,Number=1,Type=Integer,Description=\"Maximum size of the event calculated based on insert size distribution."
			"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_size_tag, &len));

	const char* ks_pval_tag = "##INFO=<ID=KS_PVAL,Number=1,Type=Float,Description=\"p-value of the KS test. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ks_pval_tag, &len));

	const char* est_size_tag = "##INFO=<ID=EST_SIZE,Number=1,Type=Integer,Description=\"Estimated size of the imprecise event. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, est_size_tag, &len));

	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_size_tag, &len));
	const char* remap_lb_tag = "##INFO=<ID=REMAP_LB,Number=1,Type=Integer,Description=\"Minimum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_lb_tag, &len));

	const char* remap_ub_tag = "##INFO=<ID=REMAP_UB,Number=1,Type=Integer,Description=\"Maximum coordinate according to the mates of the clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_ub_tag, &len));

//	const char* depths_tag = "##INFO=<ID=DEPTHS,Number=4,Type=Integer,Description=\"Depths of, respectively, the region flanking the indel to the left,"
//			"the left portion of the indel, the right portion of the indel, the region flanking the indel to the right. Numbers 2 and 3 will be identical for short indels.\">";
//	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, depths_tag, &len));

	const char* median_depths_tag = "##INFO=<ID=MEDIAN_DEPTHS,Number=4,Type=Integer,Description=\"Depths of, respectively, the region flanking the indel to the left,"
			"the left portion of the indel, the right portion of the indel, the region flanking the indel to the right. Numbers 2 and 3 will be identical for short indels.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, median_depths_tag, &len));

	const char* cluster_depths_tag = "##INFO=<ID=CLUSTER_DEPTHS,Number=2,Type=Integer,Description=\"Depths of the left and right cluster regions.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cluster_depths_tag, &len));

	const char* disc_pairs_tag = "##INFO=<ID=DISC_PAIRS,Number=1,Type=Integer,Description=\"Discordant pairs supporting the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_tag, &len));

	const char* disc_pairs_surr_tag = "##INFO=<ID=DISC_PAIRS_SURROUNDING,Number=2,Type=Integer,Description=\"Discordant pairs around the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_surr_tag, &len));

	const char* conc_pairs_tag = "##INFO=<ID=CONC_PAIRS,Number=1,Type=Integer,Description=\"Concordant pairs supporting the absence of a SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, conc_pairs_tag, &len));

	const char* clipped_reads_tag = "##INFO=<ID=CLIPPED_READS,Number=2,Type=Integer,Description=\"Reads supporting the right and the left breakpoints, respectively.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, clipped_reads_tag, &len));

	const char* max_mapq_tag = "##INFO=<ID=MAX_MAPQ,Number=2,Type=Integer,Description=\"Maximum MAPQ of a supporting clipped read.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_mapq_tag, &len));

	const char* weak_support_reads_tag = "##INFO=<ID=WEAK_SUPPORTING_1SR_READS,Number=1,Type=Integer,Description=\"Reads supporting the 1SR remapping.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_support_reads_tag, &len));

	const char* strong_support_reads_tag = "##INFO=<ID=STRONG_SUPPORTING_1SR_READS,Number=1,Type=Integer,Description=\"Reads strongly supporting the 1SR remapping.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, strong_support_reads_tag, &len));

	const char* perfect_support_reads_tag = "##INFO=<ID=PERFECT_SUPPORTING_1SR_READS,Number=1,Type=Integer,Description=\"Reads perfectly supporting the 1SR remapping.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, perfect_support_reads_tag, &len));

	const char* mmrate_pairs_tag = "##INFO=<ID=MM_RATE,Number=1,Type=Float,Description=\"Mismatch rate in consensus overlap.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, mmrate_pairs_tag, &len));

	const char* fullj_score_tag = "##INFO=<ID=FULL_JUNCTION_SCORE,Number=1,Type=Integer,Description=\"Full junction score.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fullj_score_tag, &len));

	const char* splitj_score_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE,Number=1,Type=Integer,Description=\"Split junction score.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score_tag, &len));

	const char* source_tag = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source algorithm of the indel.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, source_tag, &len));

	const char* svinsseq_tag = "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinsseq_tag, &len));

	const char* imprecise_tag = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"The reported boundaries are not precise.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, imprecise_tag, &len));

	// TODO: temporary
	const char* overlap_tag = "##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Overlap.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, overlap_tag, &len));

	const char* mismatches_tag = "##INFO=<ID=MISMATCHES,Number=1,Type=Integer,Description=\"Mismatches.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, mismatches_tag, &len));

	const char* og_range_tag = "##INFO=<ID=ORIGINAL_RANGE,Number=1,Type=String,Description=\"Unadjusted imprecise range predicted by discordant pairs. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, og_range_tag, &len));

//	// TODO: this is "debugging" information, remove when done
//	const char* extrainfo_tag = "##INFO=<ID=EXTRA_INFO,Number=.,Type=String,Description=\"Extra information.\">";
//	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, extrainfo_tag, &len));

	// add FORMAT tags
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, gt_tag, &len));

	const char* ft_tag = "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter. PASS indicates a reliable call. Any other value means the call is not reliable.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ft_tag, &len));

	// add ALT
	const char* del_alt_tag = "##ALT=<ID=DEL,Description=\"Deletion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, del_alt_tag, &len));

	const char* dup_alt_tag = "##ALT=<ID=DUP,Description=\"Tandem Duplication\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dup_alt_tag, &len));

	std::string cmd_tag = "##SurVIndel2Command=" + command;
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cmd_tag.c_str(), &len));

	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::string version_tag = "##SurVIndel2Version=" + config.version + "; Date=" + std::ctime(&now_time);
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, version_tag.c_str(), &len));

	std::stringstream called_by_ss;
	called_by_ss << "##calledBy=SurVIndel2 " << config.version << "; ";
	called_by_ss << "seed: " << config.seed << "; ";
	called_by_ss << "min_sv_size: " << config.min_sv_size << "; ";
	called_by_ss << "min_clip_len: " << config.min_clip_len << "; ";
	called_by_ss << "max_seq_error: " << config.max_seq_error << "; ";
	called_by_ss << "max_clipped_pos_dist: " << config.max_clipped_pos_dist << "; ";
	called_by_ss << "min_size_for_depth_filtering: " << config.min_size_for_depth_filtering << "; ";
	called_by_ss << "sampling-regions: " << (config.sampling_regions.empty() ? "no" : config.sampling_regions) << "; ";
	std::string called_by = called_by_ss.str();
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, called_by.c_str(), &len));

	// add samples
	bcf_hdr_add_sample(header, sample_name.c_str());

	return header;
}

void del2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, char* chr_seq, std::string& contig_name, deletion_t* del, std::vector<std::string>& filters) {
	bcf_clear(bcf_entry);
	// populate basic info
	bcf_entry->rid = bcf_hdr_name2id(hdr, contig_name.c_str());
	bcf_entry->pos = del->start;
	bcf_update_id(hdr, bcf_entry, del->id.c_str());
	std::string alleles = std::string(1, chr_seq[del->start]) + ",<DEL>";
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	// add filters
	for (std::string& filter : filters) {
		int filter_id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter.c_str());
		bcf_add_filter(hdr, bcf_entry, filter_id);
	}

	// add info
	int int_conv; // current htslib does not support writing int64 yet

	int_conv = del->end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);
	int_conv = -(del->end - del->start - del->ins_seq.length());
	bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
	bcf_update_info_string(hdr, bcf_entry, "SVTYPE", "DEL");
	if (del->max_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
		bcf_update_info_int32(hdr, bcf_entry, "MAX_SIZE", &(del->max_conf_size), 1);
	}
	if (del->remap_boundary_lower != deletion_t::REMAP_LB_NOT_COMPUTED) {
		int_conv = del->remap_boundary_lower;
		bcf_update_info_int32(hdr, bcf_entry, "REMAP_LB", &int_conv, 1);
	}
	if (del->remap_boundary_upper != deletion_t::REMAP_UB_NOT_COMPUTED) {
		int_conv = del->remap_boundary_upper;
		bcf_update_info_int32(hdr, bcf_entry, "REMAP_UB", &int_conv, 1);
	}
	int depths[] = {del->left_flanking_cov, del->indel_left_cov, del->indel_right_cov, del->right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "DEPTHS", depths, 4);
	int median_depths[] = {del->med_left_flanking_cov, del->med_indel_left_cov, del->med_indel_right_cov, del->med_right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "MEDIAN_DEPTHS", median_depths, 4);
	int cluster_depths[] = {del->med_left_cluster_cov, del->med_right_cluster_cov};
	bcf_update_info_int32(hdr, bcf_entry, "CLUSTER_DEPTHS", cluster_depths, 2);
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS", &del->disc_pairs, 1);
	int disc_pairs_surr[] = {del->l_cluster_region_disc_pairs, del->r_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
	bcf_update_info_int32(hdr, bcf_entry, "CONC_PAIRS", &del->conc_pairs, 1);
	int clipped_reads[] = {del->rc_consensus ? (int) del->rc_consensus->supp_clipped_reads : 0, del->lc_consensus ? (int) del->lc_consensus->supp_clipped_reads : 0};
	bcf_update_info_int32(hdr, bcf_entry, "CLIPPED_READS", clipped_reads, 2);
	int max_mapq[] = {del->rc_consensus ? (int) del->rc_consensus->max_mapq : 0, del->lc_consensus ? (int) del->lc_consensus->max_mapq : 0};
	bcf_update_info_int32(hdr, bcf_entry, "MAX_MAPQ", max_mapq, 2);
	if (del->sr_remap_info != NULL) {
		bcf_update_info_int32(hdr, bcf_entry, "WEAK_SUPPORTING_1SR_READS", &del->sr_remap_info->weak_supporting_reads, 1);
		bcf_update_info_int32(hdr, bcf_entry, "STRONG_SUPPORTING_1SR_READS", &del->sr_remap_info->strong_supporting_reads, 1);
		bcf_update_info_int32(hdr, bcf_entry, "PERFECT_SUPPORTING_1SR_READS", &del->sr_remap_info->perfect_supporting_reads, 1);
	}
	bcf_update_info_int32(hdr, bcf_entry, "FULL_JUNCTION_SCORE", &del->full_junction_score, 1);
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE", &del->split_junction_score, 1);
	float mm_rate = del->mm_rate;
	bcf_update_info_float(hdr, bcf_entry, "MM_RATE", &mm_rate, 1);
	bcf_update_info_string(hdr, bcf_entry, "SOURCE", del->source.c_str());
	if (!del->ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", del->ins_seq.c_str());
	}

	bcf_update_info_string(hdr, bcf_entry, "EXTRA_INFO", del->extra_info.c_str());
	bcf_update_info_flag(hdr, bcf_entry, "IMPRECISE", "", del->imprecise());

	// add GT info
	int gt[1];
	gt[0] = bcf_gt_unphased(1);
	bcf_update_genotypes(hdr, bcf_entry, gt, 1);

	const char* ft_val = (filters[0] == "PASS") ? "PASS" : "FAIL";
	bcf_update_format_string(hdr, bcf_entry, "FT", &ft_val, 1);
}

void dup2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, char* chr_seq, std::string& contig_name, duplication_t* dup, std::vector<std::string>& filters) {
	bcf_clear(bcf_entry);
	// populate basic info
	bcf_entry->rid = bcf_hdr_name2id(hdr, contig_name.c_str());
	bcf_entry->pos = dup->start;
	bcf_update_id(hdr, bcf_entry, dup->id.c_str());
	std::string alleles = std::string(1, chr_seq[dup->start]) + ",<DUP>";
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	// add filters
	for (std::string& filter : filters) {
		int filter_id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter.c_str());
		bcf_add_filter(hdr, bcf_entry, filter_id);
	}

	// add info
	int int_conv; // current htslib does not support writing int64 yet

	int_conv = dup->end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);
	int_conv = dup->end - dup->start;
	bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
	bcf_update_info_string(hdr, bcf_entry, "SVTYPE", "DUP");
	int depths[] = {dup->left_flanking_cov, dup->indel_left_cov, dup->indel_right_cov, dup->right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "DEPTHS", depths, 4);
	int median_depths[] = {dup->med_left_flanking_cov, dup->med_indel_left_cov, dup->med_indel_right_cov, dup->med_right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "MEDIAN_DEPTHS", median_depths, 4);
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS", &dup->disc_pairs, 1);
	int disc_pairs_surr[] = {dup->rc_cluster_region_disc_pairs, dup->lc_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
	int clipped_reads[] = {dup->lc_consensus ? (int) dup->lc_consensus->supp_clipped_reads : 0, dup->rc_consensus ? (int) dup->rc_consensus->supp_clipped_reads : 0};
	bcf_update_info_int32(hdr, bcf_entry, "CLIPPED_READS", clipped_reads, 2);
	int max_mapq[] = {dup->lc_consensus ? (int) dup->lc_consensus->max_mapq : 0, dup->rc_consensus ? (int) dup->rc_consensus->max_mapq : 0};
	bcf_update_info_int32(hdr, bcf_entry, "MAX_MAPQ", max_mapq, 2);
	bcf_update_info_int32(hdr, bcf_entry, "FULL_JUNCTION_SCORE", &dup->full_junction_score, 1);
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE", &dup->split_junction_score, 1);
	bcf_update_info_string(hdr, bcf_entry, "SOURCE", dup->source.c_str());

	bcf_update_info_string(hdr, bcf_entry, "EXTRA_INFO", dup->extra_info.c_str());

	// add GT info
	int gt[1];
	gt[0] = bcf_gt_unphased(1);
	bcf_update_genotypes(hdr, bcf_entry, gt, 1);

	const char* ft_val = (filters[0] == "PASS") ? "PASS" : "FAIL";
	bcf_update_format_string(hdr, bcf_entry, "FT", &ft_val, 1);
}

template<typename T>
inline T max(T a, T b, T c) { return std::max(std::max(a,b), c); }

int score(char a, char b, int match_score, int mismatch_penalty) {
	return (a == b || a == 'N' || b == 'N') ? match_score : mismatch_penalty;
}
int* smith_waterman_gotoh(const char* ref, int ref_len, const char* read, int read_len, int match_score, int mismatch_penalty, int gap_open, int gap_extend) {
	const int INF = 1000000;

	int** dab = new int*[ref_len+1];
	int** dag = new int*[ref_len+1];
	int** dgb = new int*[ref_len+1];
	for (int i = 0; i <= ref_len; i++) {
		dab[i] = new int[read_len+1];
		dag[i] = new int[read_len+1];
		dgb[i] = new int[read_len+1];
		std::fill(dab[i], dab[i]+read_len+1, 0);
		std::fill(dag[i], dag[i]+read_len+1, 0);
		std::fill(dgb[i], dgb[i]+read_len+1, 0);
	}

	dab[0][0] = dag[0][0] = dgb[0][0] = 0;
	for (int i = 1; i <= ref_len; i++) {
		dab[i][0] = -INF;
		dag[i][0] = -INF;
		dgb[i][0] = gap_open + (i-1)*gap_extend;
	}
	for (int i = 1; i <= read_len; i++) {
		dab[0][i] = -INF;
		dag[0][i] = gap_open + (i-1)*gap_extend;
		dgb[0][i] = -INF;
	}

	for (int i = 1; i <= ref_len; i++) {
		for (int j = 1; j <= read_len; j++) {
			dab[i][j] = score(ref[i-1], read[j-1], match_score, mismatch_penalty) + max(dab[i-1][j-1], dag[i-1][j-1], dgb[i-1][j-1], 0);
			dag[i][j] = max(gap_open + dab[i][j-1], gap_extend + dag[i][j-1], gap_open + dgb[i][j-1]);
			dgb[i][j] = max(gap_open + dab[i-1][j], gap_open + dag[i-1][j], gap_extend + dgb[i-1][j]);
		}
	}

	int* prefix_scores = new int[read_len];
	std::fill(prefix_scores, prefix_scores+read_len, 0);
	for (int i = 1; i <= ref_len; i++) {
		for (int j = 1; j <= read_len; j++) {
			prefix_scores[j-1] = std::max(prefix_scores[j-1], dab[i][j]);
		}
	}

	for (int i = 0; i <= ref_len; i++) {
		delete[] dab[i];
		delete[] dag[i];
		delete[] dgb[i];
	}
	delete[] dab;
	delete[] dag;
	delete[] dgb;

	for (int i = 1; i < read_len; i++) {
		prefix_scores[i] = std::max(prefix_scores[i], prefix_scores[i-1]);
	}

	return prefix_scores;
}

void remove_marked_consensuses(std::vector<consensus_t*>& consensuses, std::vector<bool>& used) {
	for (int i = 0; i < consensuses.size(); i++) if (used[i]) consensuses[i] = NULL;
	consensuses.erase(std::remove(consensuses.begin(), consensuses.end(), (consensus_t*) NULL), consensuses.end());
}

inline bool has_Ns(char* ref, hts_pos_t start, hts_pos_t len) {
	return memchr(ref+start, 'N', len) != NULL;
}

bool is_left_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 0) {
	return cigar_int_to_op(aln.cigar[0]) == 'S' && cigar_int_to_len(aln.cigar[0]) >= min_clip_len;
}
bool is_right_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 0) {
	return cigar_int_to_op(aln.cigar[aln.cigar.size()-1]) == 'S' && cigar_int_to_len(aln.cigar[aln.cigar.size()-1]) >= min_clip_len;
}
bool is_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 0) {
	return is_left_clipped(aln, min_clip_len) || is_right_clipped(aln, min_clip_len);
}

std::string get_sv_type(bcf_hdr_t* hdr, bcf1_t* sv) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
        throw std::runtime_error("Failed to determine SVTYPE for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

int get_sv_end(bcf_hdr_t* hdr, bcf1_t* sv) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        delete[] data;
        return end-1; // return 0-based
    }

    bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
    if (size > 0) {
        int svlen = data[0];
        delete[] data;
        return sv->pos + abs(svlen);
    }

    throw std::runtime_error("SV " + std::string(sv->d.id) + "has no END or SVLEN annotation.");
}

std::string get_sv_info_str(bcf_hdr_t* hdr, bcf1_t* sv, std::string info) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, info.c_str(), &data, &len) < 0) {
        throw std::runtime_error("Failed to fetch " + info + " for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

#endif //SURVINDEL2_UTILS_H
