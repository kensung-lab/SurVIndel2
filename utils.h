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

#include "simd_macros.h"


struct config_t {

    int threads, seed;
    int min_is, max_is; // find a way to move this to stats_t
    int min_sv_size;
    int match_score;
    int read_len; // this is not exactly "config", but it is more convenient to place it here
    int min_clip_len;
    double max_seq_error;
    int max_clipped_pos_dist;
    int min_size_for_depth_filtering;
    int min_diff_hsr;
    std::string sampling_regions, version;
    bool log;

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
        min_is = std::stoi(config_params["min_is"]);
        max_is = std::stoi(config_params["max_is"]);
        min_sv_size = std::stoi(config_params["min_sv_size"]);
        match_score = std::stoi(config_params["match_score"]);
        min_clip_len = std::stoi(config_params["min_clip_len"]);
        read_len = std::stod(config_params["read_len"]);
        max_seq_error = std::stod(config_params["max_seq_error"]);
        max_clipped_pos_dist = std::stoi(config_params["max_clipped_pos_dist"]);
        min_size_for_depth_filtering = std::stoi(config_params["min_size_for_depth_filtering"]);
        min_diff_hsr = std::stoi(config_params["min_diff_hsr"]);
        sampling_regions = config_params["sampling_regions"];
        version = config_params["version"];
        log = std::stoi(config_params["log"]);
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


struct consensus_t {
    bool left_clipped;
    int contig_id;
    hts_pos_t breakpoint, start, end; // we follow the vcf conventions, i.e. this is the base "before" the breakpoint
    int clip_len, lowq_clip_portion;
    std::string consensus;
    size_t supp_clipped_reads;
    uint8_t max_mapq;
    hts_pos_t remap_boundary;
    int left_ext_reads = 0, right_ext_reads = 0, hq_left_ext_reads = 0, hq_right_ext_reads = 0;
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

    hts_pos_t left_ext_target_start(config_t& config) {
    	if (!left_clipped) {
    		return start - config.max_is + config.read_len;
    	} else {
    		if (remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint - config.max_is - 2*consensus.length();
			} else {
				return remap_boundary;
			}
    	}
    }
    hts_pos_t left_ext_target_end(config_t& config) {
		if (!left_clipped) {
			return start;
		} else {
			if (remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint + config.max_is + 2*consensus.length();
			} else {
				return remap_boundary + config.max_is;
			}
		}
	}

    hts_pos_t right_ext_target_start(config_t& config) {
		if (!left_clipped) {
			if (remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
				return breakpoint - config.max_is - 2*consensus.length();
			} else {
				return remap_boundary - config.max_is;
			}
		} else {
			return end;
		}
	}
    hts_pos_t right_ext_target_end(config_t& config) {
    	if (!left_clipped) {
			if (remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
				return breakpoint + config.max_is + 2*consensus.length();
			} else {
				return remap_boundary;
			}
		} else {
			return end + config.max_is - config.read_len;
		}
    }

    int anchor_len() { return consensus.length() - clip_len; }
};

struct indel_t {
    std::string id;
    hts_pos_t start, end;
    hts_pos_t rc_anchor_start, lc_anchor_end; // start of the rc anchor and lc anchor end
    int disc_pairs = 0, disc_pairs_trusted = 0, disc_pairs_high_mapq = 0, disc_pairs_maxmapq = 0;
    consensus_t* lc_consensus,* rc_consensus;
    double mm_rate = 0.0;
    std::string source;
    int med_left_flanking_cov = 0, med_indel_left_cov = 0, med_indel_right_cov = 0, med_right_flanking_cov = 0;
    int med_left_cluster_cov = 0, med_right_cluster_cov = 0;
    int full_junction_score = 0, lh_best1_junction_score = 0, rh_best1_junction_score = 0,
    	lh_best2_junction_score = 0, rh_best2_junction_score = 0;
    int lh_junction_size = 0, rh_junction_size = 0;
    std::string extra_info;
    int overlap = 0, mismatches = 0;
    bool remapped = false;
    std::string rightmost_rightfacing_seq, leftmost_leftfacing_seq;
    std::string ins_seq;

    indel_t(hts_pos_t start, hts_pos_t end, hts_pos_t rc_anchor_start, hts_pos_t lc_anchor_end, std::string source,
            consensus_t* lc_consensus, consensus_t* rc_consensus, std::string ins_seq)
    : start(start), end(end), rc_anchor_start(rc_anchor_start), lc_anchor_end(lc_anchor_end),
	  source(source), lc_consensus(lc_consensus), rc_consensus(rc_consensus), ins_seq(ins_seq) { }

    virtual hts_pos_t len() { return end-start; }
    virtual std::string indel_type() { return ""; };
    virtual std::string to_string(std::string contig_name) {
    	return contig_name + ":" + std::to_string(start) + "-" + std::to_string(end);
    }

    bool is_single_consensus() { return (lc_consensus == NULL || rc_consensus == NULL) && lc_consensus != rc_consensus; }
    bool imprecise() { return lc_consensus == NULL && rc_consensus == NULL && remapped == false; }
};

struct deletion_t : indel_t {
	static const int SIZE_NOT_COMPUTED = INT32_MAX;
	static const int REMAP_LB_NOT_COMPUTED = 0, REMAP_UB_NOT_COMPUTED = INT32_MAX;

    int lh_score, rh_score;
    int max_conf_size = SIZE_NOT_COMPUTED, estimated_size = SIZE_NOT_COMPUTED, conc_pairs = 0;
    double ks_pval = -1.0;
    hts_pos_t remap_boundary_lower = REMAP_LB_NOT_COMPUTED, remap_boundary_upper = REMAP_UB_NOT_COMPUTED;
    int l_cluster_region_disc_pairs = 0, r_cluster_region_disc_pairs = 0;

    std::string original_range;
    std::string genotype;

    deletion_t(hts_pos_t start, hts_pos_t end, hts_pos_t rc_anchor_start, hts_pos_t lc_anchor_end, consensus_t* lc_consensus,
               consensus_t* rc_consensus, int lh_score, int rh_score, std::string source, std::string ins_seq) :
            indel_t(start, end, rc_anchor_start, lc_anchor_end, source, lc_consensus, rc_consensus, ins_seq), lh_score(lh_score), rh_score(rh_score) {};

    hts_pos_t len() override { return end-start-ins_seq.length(); }

    std::string indel_type() override { return "DEL"; }
};

struct duplication_t : indel_t {
    static const int OW_NOT_COMPUTED = INT32_MAX;

    hts_pos_t original_start, original_end;
    int lc_cluster_region_disc_pairs = 0, rc_cluster_region_disc_pairs = 0;

    duplication_t(hts_pos_t start, hts_pos_t end, hts_pos_t rc_anchor_start, hts_pos_t lc_anchor_end, consensus_t* lc_consensus,
                  consensus_t* rc_consensus, std::string source, std::string ins_seq) :
                	  original_start(start), original_end(end),
					  indel_t(std::max(start, hts_pos_t(0)), end, rc_anchor_start, lc_anchor_end, source, lc_consensus, rc_consensus, ins_seq) {}
    // Normally we store/report the base BEFORE the event, as per VCF specifications.
    // We handle here the exception where the duplication starts at the first base of the contig (base before 0 is -1, but it is not valid,
    // so we still store 0)

    std::string indel_type() override { return "DUP"; }

    hts_pos_t len() override { return ins_seq.length() + (end-start); }
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

suffix_prefix_aln_t aln_suffix_prefix_perfect(std::string& s1, std::string& s2, int min_overlap = 1) {
    int best_overlap = 0, overlap = 0;

    const char* _s1 = s1.c_str(),* _s2 = s2.c_str();
    int _s1_len = s1.length(), _s2_len = s2.length();

    int max_overlap = std::min(s1.length(), s2.length());
    for (int i = max_overlap; i >= min_overlap; i--) {
    	if (strncmp(_s1+(_s1_len-i), _s2, i) == 0) {
			return suffix_prefix_aln_t(i, i, 0);
    	}
    }
    return suffix_prefix_aln_t(0, 0, 0);
}

bool is_homopolymer(const char* seq, int len) {
	int a = 0, c = 0, g = 0, t = 0;
	for (int i = 0; i < len; i++) {
		char b = std::toupper(seq[i]);
		if (b == 'A') a++;
		else if (b == 'C') c++;
		else if (b == 'G') g++;
		else if (b == 'T') t++;
	}
	return max(a, c, g, t)/double(a+c+g+t) >= 0.8;
}

bool is_homopolymer(std::string seq) {
	return is_homopolymer(seq.data(), seq.length());
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
	const char* small_flt_tag = "##FILTER=<ID=SMALL,Description=\"Event is smaller than what required by the user.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, small_flt_tag, &len));

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

	const char* failed_to_ext_flt_tag = "##FILTER=<ID=FAILED_TO_EXTEND,Description=\"No reads can extend the consensus.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, failed_to_ext_flt_tag, &len));

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

	const char* disc_pairs_trusted_tag = "##INFO=<ID=DISC_PAIRS_TRUSTED,Number=1,Type=Integer,Description=\"Discordant pairs likely to be aligned to the correct location.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_trusted_tag, &len));

	const char* disc_pairs_hmapq_tag = "##INFO=<ID=DISC_PAIRS_HIGHMAPQ,Number=1,Type=Integer,Description=\"HDiscordant pairs with high MAPQ supporting the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_hmapq_tag, &len));

	const char* disc_pairs_surr_tag = "##INFO=<ID=DISC_PAIRS_SURROUNDING,Number=2,Type=Integer,Description=\"Discordant pairs around the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_pairs_surr_tag, &len));

	const char* conc_pairs_tag = "##INFO=<ID=CONC_PAIRS,Number=1,Type=Integer,Description=\"Concordant pairs supporting the absence of a SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, conc_pairs_tag, &len));

	const char* clipped_reads_tag = "##INFO=<ID=CLIPPED_READS,Number=2,Type=Integer,Description=\"Reads supporting the right and the left breakpoints, respectively.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, clipped_reads_tag, &len));

	const char* max_mapq_tag = "##INFO=<ID=MAX_MAPQ,Number=2,Type=Integer,Description=\"Maximum MAPQ of clipped reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, max_mapq_tag, &len));

	const char* dp_max_mapq_tag = "##INFO=<ID=DISC_PAIRS_MAXMAPQ,Number=1,Type=Integer,Description=\"Maximum MAPQ of supporting discordant pairs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dp_max_mapq_tag, &len));

	const char* ext_1sr_reads_tag = "##INFO=<ID=EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads extending a 1SR consensus to the left and to the right.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ext_1sr_reads_tag, &len));

	const char* hq_ext_1sr_reads_tag = "##INFO=<ID=HQ_EXT_1SR_READS,Number=2,Type=Integer,Description=\"Reads with high MAPQ extending a 1SR consensus to the left and to the right.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, hq_ext_1sr_reads_tag, &len));

	const char* mmrate_pairs_tag = "##INFO=<ID=MM_RATE,Number=1,Type=Float,Description=\"Mismatch rate in consensus overlap.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, mmrate_pairs_tag, &len));

	const char* fullj_score_tag = "##INFO=<ID=FULL_JUNCTION_SCORE,Number=1,Type=Integer,Description=\"Full junction score.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fullj_score_tag, &len));

	const char* splitj_score_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE,Number=2,Type=Integer,Description=\"Score of the best alignment of the left-half and right-half of the junction.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score_tag, &len));

	const char* splitj_score2_tag = "##INFO=<ID=SPLIT_JUNCTION_SCORE2,Number=2,Type=Integer,Description=\"Score of the second best alignment of the left-half and right-half of the junction.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_score2_tag, &len));

	const char* splitj_size_tag = "##INFO=<ID=SPLIT_JUNCTION_SIZE,Number=2,Type=Integer,Description=\"Size of the the left-half and right-half of the junction.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, splitj_size_tag, &len));

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

	const char* sr_consensus_seq_tag = "##INFO=<ID=SR_CONSENSUS_SEQ,Number=1,Type=String,Description=\".\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, sr_consensus_seq_tag, &len));

//	// TODO: this is "debugging" information, remove when done
	const char* extrainfo_tag = "##INFO=<ID=EXTRA_INFO,Number=.,Type=String,Description=\"Extra information.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, extrainfo_tag, &len));

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
	int median_depths[] = {del->med_left_flanking_cov, del->med_indel_left_cov, del->med_indel_right_cov, del->med_right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "MEDIAN_DEPTHS", median_depths, 4);
	int cluster_depths[] = {del->med_left_cluster_cov, del->med_right_cluster_cov};
	bcf_update_info_int32(hdr, bcf_entry, "CLUSTER_DEPTHS", cluster_depths, 2);
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS", &del->disc_pairs, 1);
	if (del->disc_pairs > 0) {
		bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_MAXMAPQ", &del->disc_pairs_maxmapq, 1);
	}
	if (del->source == "DP") {
		bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_TRUSTED", &del->disc_pairs_trusted, 1);
		bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_HIGHMAPQ", &del->disc_pairs_high_mapq, 1);
	}
	int disc_pairs_surr[] = {del->l_cluster_region_disc_pairs, del->r_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
	bcf_update_info_int32(hdr, bcf_entry, "CONC_PAIRS", &del->conc_pairs, 1);
	int clipped_reads[] = {del->rc_consensus ? (int) del->rc_consensus->supp_clipped_reads : 0, del->lc_consensus ? (int) del->lc_consensus->supp_clipped_reads : 0};
	bcf_update_info_int32(hdr, bcf_entry, "CLIPPED_READS", clipped_reads, 2);
	int max_mapq[] = {del->rc_consensus ? (int) del->rc_consensus->max_mapq : 0, del->lc_consensus ? (int) del->lc_consensus->max_mapq : 0};
	bcf_update_info_int32(hdr, bcf_entry, "MAX_MAPQ", max_mapq, 2);
	if (del->is_single_consensus()) {
		if (del->lc_consensus) {
			int ext_1sr_reads[] = { del->lc_consensus->left_ext_reads, del->lc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { del->lc_consensus->hq_left_ext_reads, del->lc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
			bcf_update_info_string(hdr, bcf_entry, "SR_CONSENSUS_SEQ", del->lc_consensus->consensus.c_str());
		} else if (del->rc_consensus) {
			int ext_1sr_reads[] = { del->rc_consensus->left_ext_reads, del->rc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { del->rc_consensus->hq_left_ext_reads, del->rc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
			bcf_update_info_string(hdr, bcf_entry, "SR_CONSENSUS_SEQ", del->rc_consensus->consensus.c_str());
		}
	}
	bcf_update_info_int32(hdr, bcf_entry, "FULL_JUNCTION_SCORE", &del->full_junction_score, 1);
	int split_junction_score[] = {del->lh_best1_junction_score, del->rh_best1_junction_score};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE", split_junction_score, 2);
	int split_junction_score2[] = {del->lh_best2_junction_score, del->rh_best2_junction_score};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE2", split_junction_score2, 2);
	int split_junction_size[] = {del->lh_junction_size, del->rh_junction_size};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SIZE", split_junction_size, 2);
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
	std::string alleles = std::string(1, chr_seq[dup->start]);
	if (dup->start == dup->end) {
		alleles += ",<INS>";
	} else {
		alleles += ",<DUP>";
	}
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
	int_conv = dup->len();
	bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
	if (dup->start == dup->end) {
		bcf_update_info_string(hdr, bcf_entry, "SVTYPE", "INS");
	} else {
		bcf_update_info_string(hdr, bcf_entry, "SVTYPE", "DUP");
	}
	int median_depths[] = {dup->med_left_flanking_cov, dup->med_indel_left_cov, dup->med_indel_right_cov, dup->med_right_flanking_cov};
	bcf_update_info_int32(hdr, bcf_entry, "MEDIAN_DEPTHS", median_depths, 4);
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS", &dup->disc_pairs, 1);
	int disc_pairs_surr[] = {dup->rc_cluster_region_disc_pairs, dup->lc_cluster_region_disc_pairs};
	bcf_update_info_int32(hdr, bcf_entry, "DISC_PAIRS_SURROUNDING", disc_pairs_surr, 2);
	int clipped_reads[] = {dup->lc_consensus ? (int) dup->lc_consensus->supp_clipped_reads : 0, dup->rc_consensus ? (int) dup->rc_consensus->supp_clipped_reads : 0};
	bcf_update_info_int32(hdr, bcf_entry, "CLIPPED_READS", clipped_reads, 2);
	int max_mapq[] = {dup->lc_consensus ? (int) dup->lc_consensus->max_mapq : 0, dup->rc_consensus ? (int) dup->rc_consensus->max_mapq : 0};
	bcf_update_info_int32(hdr, bcf_entry, "MAX_MAPQ", max_mapq, 2);
	if (dup->is_single_consensus()) {
		if (dup->lc_consensus) {
			int ext_1sr_reads[] = { dup->lc_consensus->left_ext_reads, dup->lc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { dup->lc_consensus->hq_left_ext_reads, dup->lc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
			bcf_update_info_string(hdr, bcf_entry, "SR_CONSENSUS_SEQ", dup->lc_consensus->consensus.c_str());
		} else if (dup->rc_consensus) {
			int ext_1sr_reads[] = { dup->rc_consensus->left_ext_reads, dup->rc_consensus->right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "EXT_1SR_READS", ext_1sr_reads, 2);
			int hq_ext_1sr_reads[] = { dup->rc_consensus->hq_left_ext_reads, dup->rc_consensus->hq_right_ext_reads };
			bcf_update_info_int32(hdr, bcf_entry, "HQ_EXT_1SR_READS", hq_ext_1sr_reads, 2);
			bcf_update_info_string(hdr, bcf_entry, "SR_CONSENSUS_SEQ", dup->rc_consensus->consensus.c_str());
		}
	}
	bcf_update_info_int32(hdr, bcf_entry, "FULL_JUNCTION_SCORE", &dup->full_junction_score, 1);
	int split_junction_score[] = {dup->lh_best1_junction_score, dup->rh_best1_junction_score};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE", split_junction_score, 2);
	int split_junction_score2[] = {dup->lh_best2_junction_score, dup->rh_best2_junction_score};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SCORE2", split_junction_score2, 2);
	int split_junction_size[] = {dup->lh_junction_size, dup->rh_junction_size};
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_JUNCTION_SIZE", split_junction_size, 2);
	bcf_update_info_string(hdr, bcf_entry, "SOURCE", dup->source.c_str());

	if (!dup->ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", dup->ins_seq.c_str());
	}
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

size_t gcd(size_t a, size_t b) {
    return (b == 0) ? a : gcd(b, a % b);
}
size_t lcm(size_t a, size_t b) {
    return a * b / gcd(a, b);
}
int* smith_waterman_gotoh(const char* ref, int ref_len, const char* read, int read_len,
		int match_score, int mismatch_penalty, int gap_open, int gap_extend) {

	const int INF = 1000000;

	const int BYTES_PER_BLOCK = INT_PER_BLOCK * 4;

	// turn read_len+1 into a multiple of INT_PER_BLOCK
	int read_len_rounded = (read_len+1+INT_PER_BLOCK-1)/INT_PER_BLOCK*INT_PER_BLOCK;

	const char* alphabet = "NACGT";
	const int alphabet_size = strlen(alphabet);

	int* H = NULL;
	int* E = NULL;
	int* F = NULL;
	int* prefix_scores = NULL;

	int** profile = new int*[alphabet_size];
	size_t alignment = lcm(BYTES_PER_BLOCK, sizeof(void*));
	int p1 = posix_memalign(reinterpret_cast<void**>(&H), alignment, 2*read_len_rounded * sizeof(int));
	int p2 = posix_memalign(reinterpret_cast<void**>(&E), alignment, 2*read_len_rounded * sizeof(int));
	int p3 = posix_memalign(reinterpret_cast<void**>(&F), alignment, 2*read_len_rounded * sizeof(int));
	int p4 = posix_memalign(reinterpret_cast<void**>(&prefix_scores), alignment, read_len_rounded * sizeof(int));
	int p5 = 0;
	for (int i = 0; i < alphabet_size; i++) {
		p5 += posix_memalign(reinterpret_cast<void**>(&profile[i]), alignment, read_len_rounded * sizeof(int));
		profile[i][0] = 0;
		for (int j = 1; j <= read_len; j++) {
			profile[i][j] = (read[j-1] == alphabet[i]) ? match_score : mismatch_penalty;
		}
		for (int j = read_len+1; j < read_len_rounded; j++) {
			profile[i][j] = 0;
		}
	}
	if (p1 || p2 || p3 || p4 || p5) {
		std::cerr << "Error allocating aligned memory of size " << (2*read_len_rounded * sizeof(int)) << std::endl;
	}

	int* H_prev = H, *H_curr = H+read_len_rounded;
	int* E_prev = E, *E_curr = E+read_len_rounded;
	int* F_prev = F, *F_curr = F+read_len_rounded;

	std::fill(H_prev, H_prev+read_len_rounded, 0);
	std::fill(E_prev, E_prev+read_len_rounded, 0);
	std::fill(F_prev, F_prev+read_len_rounded, 0);
	H_curr[0] = F_curr[0] = E_curr[0] = 0;

	SIMD_INT gap_open_v = SET1_INT(gap_open);
	SIMD_INT gap_open_v_pos = SET1_INT(-gap_open);
	SIMD_INT gap_extend_v = SET1_INT(gap_extend);
	SIMD_INT zero_v = SET1_INT(0);

	std::fill(prefix_scores, prefix_scores+read_len, 0);
	for (int i = 1; i <= ref_len; i++) {
		for (int j = 0; j < read_len_rounded; j += INT_PER_BLOCK) {
			SIMD_INT H_up_v = LOAD_INT((SIMD_INT*)&H_prev[j]);
			SIMD_INT E_up_v = LOAD_INT((SIMD_INT*)&E_prev[j]);
			SIMD_INT F_up_v = LOAD_INT((SIMD_INT*)&F_prev[j]);
			SIMD_INT m1 = ADD_INT(gap_open_v, MAX_INT(H_up_v, F_up_v));
			SIMD_INT E_curr_v = MAX_INT(m1, ADD_INT(gap_extend_v, E_up_v));
			STORE_INT((SIMD_INT*)&E_curr[j], E_curr_v);
		}

		int* ref_profile = profile[0];
		switch (ref[i-1]) {
			case 'A': ref_profile = profile[1]; break;
			case 'C': ref_profile = profile[2]; break;
			case 'G': ref_profile = profile[3]; break;
			case 'T': ref_profile = profile[4]; break;
		}
		for (int j = 1; j <= read_len_rounded-INT_PER_BLOCK; j += INT_PER_BLOCK) {
			SIMD_INT H_diag_v = LOAD_INT((SIMD_INT*)&H_prev[j-1]);
			SIMD_INT F_diag_v = LOAD_INT((SIMD_INT*)&F_prev[j-1]);
			SIMD_INT E_diag_v = LOAD_INT((SIMD_INT*)&E_prev[j-1]);
			SIMD_INT m1 = MAX_INT(H_diag_v, F_diag_v);
			m1 = MAX_INT(m1, E_diag_v);
			SIMD_INT H_curr_v = MAX_INT(m1, zero_v);

			SIMD_INT profile_curr = LOADU_INT((SIMD_INT*)&ref_profile[j]);
			H_curr_v = ADD_INT(H_curr_v, profile_curr);

			STOREU_INT((SIMD_INT*)&H_curr[j], H_curr_v);
			
			SIMD_INT prefix_v = LOAD_INT((SIMD_INT*)&prefix_scores[j-1]);
			prefix_v = MAX_INT(prefix_v, H_curr_v);
			STORE_INT((SIMD_INT*)&prefix_scores[j-1], prefix_v);
		}
		for (int j = read_len_rounded-INT_PER_BLOCK; j < read_len_rounded; j++) {
			H_curr[j] = ref_profile[j] + max(H_prev[j-1], F_prev[j-1], E_prev[j-1], 0);
			prefix_scores[j] = std::max(prefix_scores[j], H_curr[j]);
		}

		int j = 0;
		for (; j < read_len_rounded; j += INT_PER_BLOCK) {
			SIMD_INT H_curr_v = LOAD_INT((SIMD_INT*)&H_curr[j]);
			auto cmp = CMP_INT(H_curr_v, gap_open_v_pos);
			if (cmp) break;
			STORE_INT((SIMD_INT*)&F_curr[j], zero_v);
		}
		if (j == 0) j = 1;
		for (; j <= read_len; j++) {
			F_curr[j] = std::max(gap_open + H_curr[j-1], gap_extend + F_curr[j-1]);
		}

		std::swap(H_prev, H_curr);
		std::swap(E_prev, E_curr);
		std::swap(F_prev, F_curr);
	}

	free(H);
	free(E);
	free(F);
	for (int i = 0; i < alphabet_size; i++) free(profile[i]);
	delete[] profile;

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

int get_left_clip_size(StripedSmithWaterman::Alignment& aln) {
	return cigar_int_to_op(aln.cigar[0]) == 'S' ? cigar_int_to_len(aln.cigar[0]) : 0;
}
int get_right_clip_size(StripedSmithWaterman::Alignment& aln) {
	return cigar_int_to_op(aln.cigar[aln.cigar.size()-1]) == 'S' ? cigar_int_to_len(aln.cigar[aln.cigar.size()-1]) : 0;
}

bool is_left_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return get_left_clip_size(aln) >= min_clip_len;
}
bool is_right_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return get_right_clip_size(aln) >= min_clip_len;
}
bool is_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return is_left_clipped(aln, min_clip_len) || is_right_clipped(aln, min_clip_len);
}

std::pair<int, int> find_aln_prefix_score(std::vector<uint32_t> cigar, int ref_prefix_len, int match_score, int mismatch_score,
						  int gap_open_score, int gap_extend_score) {
	int score = 0, query_prefix_len = 0;
	for (int i = 0, j = 0; i < cigar.size() && j < ref_prefix_len; i++) {
		int op = cigar_int_to_op(cigar[i]);
		int len = cigar_int_to_len(cigar[i]);
		if (op == '=') {
			len = std::min(len, ref_prefix_len-j);
			score += len*match_score;
			query_prefix_len += len;
			j += len;
		} else if (op == 'X') {
			len = std::min(len, ref_prefix_len-j);
			score += len*mismatch_score;
			query_prefix_len += len;
			j += len;
		} else if (op == 'I') {
			score += gap_open_score + (len-1)*gap_extend_score;
			query_prefix_len += len;
		} else if (op == 'D') {
			len = std::min(len, ref_prefix_len-j);
			score += gap_open_score + (len-1)*gap_extend_score;
			j += len;
		}
	}
	return {score, query_prefix_len};
}

std::pair<int, int> find_aln_suffix_score(std::vector<uint32_t> cigar, int ref_suffix_len, int match_score, int mismatch_score,
						  int gap_open_score, int gap_extend_score) {
	std::vector<uint32_t> rev_cigar(cigar.rbegin(), cigar.rend());
	return find_aln_prefix_score(rev_cigar, ref_suffix_len, match_score, mismatch_score, gap_open_score, gap_extend_score);
}

#endif //SURVINDEL2_UTILS_H
