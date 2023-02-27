#include <string>
#include <cmath>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "sam_utils.h"
#include "utils.h"
//#include "hsr.h"
#include "libs/cptl_stl.h"

config_t config;

std::string workdir, reference_fname;

std::vector<std::pair<uint64_t, uint32_t> > partial_sums;
std::vector<uint32_t> depths;
std::vector<int> as_diff_dist;
std::mutex mtx;

uint32_t* isize_counts;
std::vector<std::vector<uint32_t> > isizes_count_geq_i;

bool is_hidden_split_read(bam1_t* r) {
	if (is_left_clipped(r, 0) || is_right_clipped(r, 0)) return false;

    int mismatches = bam_aux2i(bam_aux_get(r, "NM"));

    int indels = 0;
    uint32_t* cigar = bam_get_cigar(r);
    for (uint32_t i = 0; i < r->core.n_cigar; i++) {
        char op_chr = bam_cigar_opchr(cigar[i]);
        if (op_chr == 'D' || op_chr == 'I') {
            mismatches -= bam_cigar_oplen(cigar[i]);
            indels++;
        }
    }

    return (mismatches + indels >= 3);
}

void categorize(int id, std::string bam_fname, int contig_id, std::string contig_name, std::vector<int> rnd_positions) {

    sort(rnd_positions.begin(), rnd_positions.end());

    open_samFile_t* bam_file = open_samFile(bam_fname);
    if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
        throw "Failed to read reference " + reference_fname;
    }

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    std::string prefix = workdir + "/workspace/" + std::to_string(contig_id);
    samFile* clip_writer = get_writer(prefix + "-CLIP.bam", bam_file->header);
    samFile* hsr_writer = get_writer(prefix + "-HSR.bam", bam_file->header);
    samFile* dp_writer = get_writer(prefix + "-DP.bam", bam_file->header);
    std::ofstream mateseqs_fout(prefix + ".mateseqs");

    uint64_t sum_is = 0;
    uint32_t n_is = 0;

    int curr_pos = 0;
    std::vector<uint32_t> rnd_positions_depth(rnd_positions.size());
    std::vector<int> as_diff_dist_contig;
    std::vector<uint32_t> local_isize_counts(config.max_is+1);
    std::vector<std::vector<uint32_t> > local_isize_dist_by_pos(rnd_positions.size());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;

        while (curr_pos < rnd_positions.size() && read->core.pos > rnd_positions[curr_pos]) curr_pos++;

        int rlen = read->core.l_qseq;
        if (curr_pos < rnd_positions.size() && is_samechr(read) && !is_samestr(read) && !bam_is_rev(read)
        && read->core.isize > 0 && read->core.isize <= config.max_is) {
            if (read->core.pos+rlen/2 <= rnd_positions[curr_pos] && rnd_positions[curr_pos] <= read->core.pos+read->core.isize-rlen/2) {
            	local_isize_counts[read->core.isize]++;
            	local_isize_dist_by_pos[curr_pos].push_back(read->core.isize);
                sum_is += read->core.isize;
                n_is++;
            }
        }

        // clipped read
        if ((is_left_clipped(read, config.min_clip_len)) ||
            (is_right_clipped(read, config.min_clip_len))) {
            int ok = sam_write1(clip_writer, bam_file->header, read);
            if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
        } else if (is_samechr(read) && is_hidden_split_read(read)) {
            int ok = sam_write1(hsr_writer, bam_file->header, read);
            if (ok < 0) throw "Failed to write to " + std::string(hsr_writer->fn);
        }

        if (is_samechr(read) && !is_samestr(read) && is_primary(read)) {
        	if (!bam_is_rev(read) && (read->core.isize > config.max_is || read->core.isize < 0)) {
				int ok = sam_write1(dp_writer, bam_file->header, read);
				if (ok < 0) throw "Failed to write to " + std::string(dp_writer->fn);
        	} else if (bam_is_rev(read) && (read->core.isize < -config.max_is || read->core.isize > 0)) {
        		mateseqs_fout << bam_get_qname(read) << " " << get_sequence(read) << std::endl;
        	}
        }

        if (!is_samechr(read) || read->core.qual < 20) continue;

        // add depth
        hts_pos_t endpos = bam_endpos(read);
        for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < endpos; i++) {
            if (read->core.pos <= rnd_positions[i] && bam_endpos(read) >= rnd_positions[i]) {
                rnd_positions_depth[i]++;
            }
        }

        // sample deviation compared to optimal score (which is read length)
        // but only if the read is VERY confidently aligned (primary, MAPQ = 60 and not clipped)
        // this is because we (as much as possible, and since actual SNPs/indels are *relatively* rare) want the deviation
        // in score to be due to actual sequencing errors
        if (read->core.qual == 60 && !is_left_clipped(read, 0) && !is_right_clipped(read, 0)) {
            int64_t as = get_AS_tag(read);
            as_diff_dist_contig.push_back(config.match_score*read->core.l_qseq - as);
        }
    }

    mtx.lock();
    partial_sums.push_back({sum_is, n_is});
    for (uint32_t d : rnd_positions_depth) {
        if (d > 0) depths.push_back(d);
    }
    for (int as_diff : as_diff_dist_contig) {
        as_diff_dist.push_back(as_diff);
    }

    for (int i = 0; i <= config.max_is; i++) {
    	isize_counts[i] += local_isize_counts[i];
    }

    for (int i = 0; i < rnd_positions.size(); i++) {
		std::vector<uint32_t> local_isizes_count_geq_i(config.max_is+1); // position i contains the count of isizes that are geq than i
		for (uint32_t isize : local_isize_dist_by_pos[i]) local_isizes_count_geq_i[isize]++;
		for (int j = config.max_is-1; j >= 0; j--) {
			local_isizes_count_geq_i[j] += local_isizes_count_geq_i[j+1];
		}
		if (local_isizes_count_geq_i[0] == 0) continue;
		for (int j = 0; j <= config.max_is; j++) {
			isizes_count_geq_i[j].push_back(local_isizes_count_geq_i[j]);
		}
    }

    mtx.unlock();

    sam_close(clip_writer);
    sam_close(hsr_writer);
    sam_close(dp_writer);

    close_samFile(bam_file);
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    reference_fname = argv[3];

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");

    std::string contig_name;
    std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
    int pos;
    std::unordered_map<std::string, std::vector<int> > rnd_pos_map;
    while (rnd_pos_fin >> contig_name >> pos) {
        rnd_pos_map[contig_name].push_back(pos);
    }

    isize_counts = new uint32_t[config.max_is+1];
    std::fill(isize_counts, isize_counts+config.max_is+1, 0);
    isizes_count_geq_i.resize(config.max_is+1);

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool.push(categorize, bam_fname, contig_id, contig_name, rnd_pos_map[contig_name]);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cout << s << std::endl;
        }
    }

    uint64_t sum_is = 0;
    uint32_t n_is = 0;
    for (std::pair<uint64_t, uint32_t>& p : partial_sums) {
        sum_is += p.first;
        n_is += p.second;
    }

    uint64_t lt_depth_stddev = 0, gt_depth_stddev = 0;
    int gt = 0, lt = 0;
    uint32_t avg_depth = mean(depths);
    for (uint32_t d : depths) {
        if (d <= avg_depth) {
            lt_depth_stddev += (avg_depth-d)*(avg_depth-d);
            lt++;
        } else if (d >= avg_depth) {
            gt_depth_stddev += (avg_depth-d)*(avg_depth-d);
            gt++;
        }
    }
    lt_depth_stddev = std::sqrt(lt_depth_stddev/(lt-1));
    gt_depth_stddev = std::sqrt(gt_depth_stddev/(gt-1));

    std::sort(depths.begin(), depths.end());

    std::ofstream stats_out(workdir + "/stats.txt");
    stats_out << "pop_avg_crossing_is " << sum_is / n_is << std::endl;
    stats_out << "avg_depth " << avg_depth << std::endl;
    stats_out << "lt_depth_stddev " << lt_depth_stddev << std::endl;
    stats_out << "gt_depth_stddev " << gt_depth_stddev << std::endl;
    stats_out << "min_depth " << depths[depths.size()/100] << std::endl;
    stats_out << "max_depth " << depths[depths.size()-depths.size()/100] << std::endl;

    // convert as_diff_dist into a frequency distribution
    std::sort(as_diff_dist.begin(), as_diff_dist.end());
    int max_as_diff = as_diff_dist[as_diff_dist.size()-1];
    std::vector<uint64_t> as_diff_counts(max_as_diff+1);
    for (int as_diff : as_diff_dist) {
        as_diff_counts[as_diff]++;
    }

    std::ofstream as_diff_fout(workdir + "/as_diff_dist.txt");
    for (int i = 0; i <= max_as_diff; i++) {
        as_diff_fout << i << " " << double(as_diff_counts[i])/as_diff_dist.size() << "\n";
    }
    as_diff_fout.close();

    std::ofstream crossing_isizes_dist_fout(workdir + "/crossing_isizes.txt");
    for (int i = 0; i <= config.max_is; i++) {
    	crossing_isizes_dist_fout << i << " " << isize_counts[i] << std::endl;
    }
    crossing_isizes_dist_fout.close();
    delete[] isize_counts;

    std::ofstream crossing_isizes_count_geq_i_fout(workdir + "/crossing_isizes_count_geq_i.txt");
    int mid = isizes_count_geq_i[0].size()/2;
	for (int i = 0; i <= config.max_is; i++) {
		crossing_isizes_count_geq_i_fout << i << " ";
		std::sort(isizes_count_geq_i[i].begin(), isizes_count_geq_i[i].end());
		crossing_isizes_count_geq_i_fout << isizes_count_geq_i[i][mid] << std::endl;
	}
	crossing_isizes_dist_fout.close();

	stats_out << "min_disc_pairs " << isizes_count_geq_i[0][isizes_count_geq_i[0].size()/100] << std::endl;
	stats_out.close();
}
