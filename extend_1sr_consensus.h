#ifndef EXTEND_1SR_CONSENSUS_H_
#define EXTEND_1SR_CONSENSUS_H_

#include <chrono>
#include <unordered_map>
#include <queue>
#include <stack>
#include "utils.h"

struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

struct ext_read_t {
	std::string qname;
	hts_pos_t start, end;
	uint8_t mapq;
	bool rev, same_chr;
	char* sequence = NULL;
	// int64_t* packed_seq = NULL;

	ext_read_t(bam1_t* read) : qname(bam_get_qname(read)), start(read->core.pos), end(bam_endpos(read)), mapq(read->core.qual),
			rev(bam_is_rev(read)), same_chr(is_samechr(read)), sequence(get_sequence_cstr(read)) {
		// std::string seq = get_sequence(read);
		// packed_seq = new int64_t[seq.length()/32 + 1];
		// memset(packed_seq, 0, sizeof(int64_t)*(seq.length()/32 + 1));
		// // store seq as 2-bit packed int64_t
		// for (int i = 0; i < seq.length(); i++) {
		// 	int64_t nv = 0;
		// 	if (seq[i] == 'A' || seq[i] == 'a') nv = 0;
		// 	else if (seq[i] == 'C' || seq[i] == 'c') nv = 1;
		// 	else if (seq[i] == 'G' || seq[i] == 'g') nv = 2;
		// 	else if (seq[i] == 'T' || seq[i] == 't') nv = 3;
		// 	else nv = 0;
		// 	packed_seq[i/32] = (packed_seq[i/32] << 2) | nv;
		// }
	}

	~ext_read_t() {
		if (sequence != NULL) delete[] sequence;
		// if (packed_seq != NULL) delete[] packed_seq;
	}
};

std::vector<bool> find_vertices_in_cycles(std::vector<std::vector<edge_t> >& l_adj) {

	const int INF = 1000000;

	int n = l_adj.size();
	std::vector<std::vector<int> > dist;
	for (int i = 0; i < n; i++) dist.push_back(std::vector<int>(n, INF));
	for (int i = 0; i < n; i++) {
		for (edge_t& e : l_adj[i]) dist[i][e.next] = 1;
	}
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (dist[i][j] > dist[i][k] + dist[k][j]) {
					dist[i][j] = dist[i][k] + dist[k][j];
				}
			}
		}
	}

	std::vector<bool> in_cycle(n);
	for (int i = 0; i < n; i++) {
		in_cycle[i] = (dist[i][i] != INF);
	}
	return in_cycle;
}

bool is_vertex_in_cycle(std::vector<std::vector<edge_t> >& l_adj, int i) {
	std::stack<int> s;
	for (edge_t& e : l_adj[i]) {
		s.push(e.next);
	}

	// check if node i can reach itself - i.e., it is part of a cycle
	std::vector<bool> visited(l_adj.size(), false);
	bool i_in_cycle = false;
	while (!s.empty()) {
		int curr = s.top();
		s.pop();
		if (visited[curr]) continue;
		visited[curr] = true;

		if (curr == i) return true;
		for (edge_t& e : l_adj[curr]) {
			s.push(e.next);
		}
	}
	return false;
}

std::vector<bool> find_vertices_in_cycles_fast(std::vector<std::vector<edge_t> >& l_adj) {

	int n = l_adj.size();
	std::vector<bool> in_cycle(n);
	for (int i = 0; i < n; i++) {
		in_cycle[i] = is_vertex_in_cycle(l_adj, i);
	}
	return in_cycle;
}

void build_graph_fwd(std::vector<std::string>& read_seqs, std::vector<int> starting_idxs, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, int min_overlap) {

	uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

	int n = read_seqs.size();
	std::vector<bool> visited(n);

	std::unordered_map<uint64_t, std::vector<int> > kmer_to_idx;
	for (int i = 0; i < n; i++) {
		uint64_t kmer = 0;

		std::string& seq = read_seqs[i];
		for (int j = 0; j < seq.length()-1; j++) {
			uint64_t nv = nucl_bm[seq[j]];
			kmer = ((kmer << 2) | nv);

			if (j >= 32) {
				kmer_to_idx[kmer].push_back(i);
			}
		}
	}

	std::queue<int> bfs;
	for (int i : starting_idxs) bfs.push(i);
	while (!bfs.empty()) {
		int i = bfs.front();
		bfs.pop();

		if (visited[i]) continue;
		visited[i] = true;

		std::string& s1 = read_seqs[i];
		uint64_t kmer = 0;
		for (int j = s1.length()-32; j < s1.length(); j++) {
			uint64_t nv = nucl_bm[s1[j]];
			kmer = ((kmer << 2) | nv);
		}

		for (int j : kmer_to_idx[kmer]) {
			std::string& s2 = read_seqs[j];

			if (s1 == s2) continue;

			suffix_prefix_aln_t spa = aln_suffix_prefix_perfect(s1, s2, min_overlap);
			if (spa.overlap) {
				if (!is_homopolymer(s2.c_str(), spa.overlap)) {
					out_edges[i]++;
					l_adj[i].push_back({j, spa.score, spa.overlap});
					l_adj_rev[j].push_back({i, spa.score, spa.overlap});
					bfs.push(j);
				}
			}
		}
	}
}

void build_graph_rev(std::vector<std::string>& read_seqs, std::vector<int> starting_idxs, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, int min_overlap) {

	uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

	int n = read_seqs.size();
	std::vector<bool> visited(n);

	std::unordered_map<uint64_t, std::vector<int> > kmer_to_idx;

	for (int i = 0; i < n; i++) {
		uint64_t kmer = 0;

		std::string& seq = read_seqs[i];
		for (int j = 0; j < seq.length(); j++) {
			uint64_t nv = nucl_bm[seq[j]];
			kmer = ((kmer << 2) | nv);

			if (j >= 32+1) {
				kmer_to_idx[kmer].push_back(i);
			}
		}
	}

	std::queue<int> bfs;
	for (int i : starting_idxs) bfs.push(i);
	while (!bfs.empty()) {
		int i = bfs.front();
		bfs.pop();

		if (visited[i]) continue;
		visited[i] = true;

		std::string& s1 = read_seqs[i];
		uint64_t kmer = 0;
		for (int j = 0; j < 32; j++) {
			uint64_t nv = nucl_bm[s1[j]];
			kmer = ((kmer << 2) | nv);
		}

		for (int j : kmer_to_idx[kmer]) {
			std::string& s2 = read_seqs[j];

			if (s1 == s2) continue;

			suffix_prefix_aln_t spa2 = aln_suffix_prefix_perfect(s2, s1, min_overlap);
			if (spa2.overlap) {
				bool spa2_homopolymer = is_homopolymer(s1.c_str(), spa2.overlap);
				if (!spa2_homopolymer) {
					out_edges[i]++;
					l_adj[i].push_back({j, spa2.score, spa2.overlap});
					l_adj_rev[j].push_back({i, spa2.score, spa2.overlap});
					bfs.push(j);
				}
			}
		}
	}
}

std::vector<int> find_rev_topological_order(int n, std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj_rev) {

	std::queue<int> sinks;
	for (int i = 0; i < n; i++) {
		if (!out_edges[i]) sinks.push(i);
	}

	std::vector<int> rev_topological_order;
	while (!sinks.empty()) {
		int s = sinks.front();
		sinks.pop();
		rev_topological_order.push_back(s);
		for (edge_t& e : l_adj_rev[s]) {
			out_edges[e.next]--;
			if (out_edges[e.next] == 0) sinks.push(e.next);
		}
	}
	return rev_topological_order;
}

void get_extension_read_seqs(IntervalTree<ext_read_t*>& candidate_reads_itree, std::vector<std::string>& read_seqs, std::vector<int>& read_mapqs,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq,
		hts_pos_t target_start, hts_pos_t target_end, hts_pos_t contig_len, config_t config, int max_reads = INT32_MAX) {

	hts_pos_t fwd_mates_start = std::max(hts_pos_t(1), target_start-config.max_is+config.read_len);
	hts_pos_t fwd_mates_end = std::max(hts_pos_t(1), target_end-config.min_is);

	hts_pos_t rev_mates_start = std::min(target_start+config.min_is, contig_len);
	hts_pos_t rev_mates_end = std::min(target_end+config.max_is-config.read_len, contig_len);

	std::vector<Interval<ext_read_t*> > target_reads = candidate_reads_itree.findOverlapping(
			std::min(fwd_mates_start, rev_mates_start)-10, std::max(fwd_mates_end, rev_mates_end)+10);
	for (Interval<ext_read_t*> i_read : target_reads) {
		ext_read_t* ext_read = i_read.value;
		if (!ext_read->rev && ext_read->mapq >= config.high_confidence_mapq &&
			!ext_read->same_chr && ext_read->start >= fwd_mates_start && ext_read->start <= fwd_mates_end) {
			std::pair<std::string, int> mateseq_w_mapq = mateseqs_w_mapq[ext_read->qname];
			rc(mateseq_w_mapq.first);
			read_seqs.push_back(mateseq_w_mapq.first);
			read_mapqs.push_back(mateseq_w_mapq.second);
		} else if (ext_read->rev && ext_read->mapq >= config.high_confidence_mapq &&
				!ext_read->same_chr && ext_read->end >= rev_mates_start && ext_read->end <= rev_mates_end) {
			std::pair<std::string, int> mateseq_w_mapq = mateseqs_w_mapq[ext_read->qname];
			read_seqs.push_back(mateseq_w_mapq.first);
			read_mapqs.push_back(mateseq_w_mapq.second);
		}

		if (ext_read->end > target_start && ext_read->start < target_end) {
			read_seqs.push_back(ext_read->sequence);
			read_mapqs.push_back(ext_read->mapq);
		}

		if (read_seqs.size() > max_reads) return;
	}

	// TODO : results very slightly change when sorting, investigate
//	std::vector<std::pair<std::string, int>> temp;
//	for (int i = 1; i < read_seqs.size(); i++) {
//		temp.push_back({read_seqs[i], read_mapqs[i]});
//	}
//	std::sort(temp.begin(), temp.end());
//
//	read_seqs = std::vector<std::string>(1, read_seqs[0]);
//	read_mapqs = std::vector<int>(1, read_mapqs[0]);
//	for (int i = 0; i < temp.size(); i++) {
//		read_seqs.push_back(temp[i].first);
//		read_mapqs.push_back(temp[i].second);
//	}
}
std::vector<ext_read_t*> get_extension_reads(std::string contig_name, std::vector<hts_pair_pos_t>& target_ivals, hts_pos_t contig_len,
		std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq, config_t config, open_samFile_t* bam_file) {

	std::vector<char*> regions;
	for (hts_pair_pos_t target_ival : target_ivals) {
		hts_pos_t fwd_mates_start = std::max(hts_pos_t(1), target_ival.beg-config.max_is+config.read_len);
		hts_pos_t rev_mates_end = std::min(target_ival.end+config.max_is-config.read_len, contig_len);

		std::stringstream ss;
		ss << contig_name << ":" << fwd_mates_start << "-" << rev_mates_end;

		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		regions.push_back(region);
	}

	std::vector<ext_read_t*> reads;

	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());

	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		reads.push_back(new ext_read_t(read));
	}

	bam_destroy1(read);
	hts_itr_destroy(iter);

	for (char* region : regions) delete[] region;

	return reads;
}

void break_cycles(std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev) {
	int n = l_adj.size();
	std::vector<bool> in_cycle = find_vertices_in_cycles_fast(l_adj);
	for (int i = 0; i < n; i++) {
		if (in_cycle[i]) {
			l_adj[i].clear();
			out_edges[i] = 0;
		}
		l_adj_rev[i].erase(
			std::remove_if(l_adj_rev[i].begin(), l_adj_rev[i].end(), [&in_cycle](edge_t& e) {return in_cycle[e.next];}),
			l_adj_rev[i].end()
		);
	}
}

void extend_consensus_to_right(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree,
		hts_pos_t target_start, hts_pos_t target_end, std::string contig_name, hts_pos_t contig_len,
		config_t config, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	std::vector<std::string> read_seqs;
	std::vector<int> read_mapqs;
	read_seqs.push_back(consensus->consensus);
	read_mapqs.push_back(config.high_confidence_mapq);

	get_extension_read_seqs(candidate_reads_itree, read_seqs, read_mapqs, mateseqs_w_mapq, target_start, target_end, contig_len, config, 5000);

	if (read_seqs.size() > 5000) return;

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	build_graph_fwd(read_seqs, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, config.read_len/2);

	bool cycle_contains_0 = is_vertex_in_cycle(l_adj, 0);
	if (cycle_contains_0) {
		// remove all edges pointing to 0
		for (std::vector<edge_t>& l_adj_i : l_adj) {
			l_adj_i.erase(std::remove_if(l_adj_i.begin(), l_adj_i.end(), [](edge_t& e) {return e.next == 0;}), l_adj_i.end());
		}
		l_adj_rev[0].clear();
	}

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	bool cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) {
		break_cycles(out_edges, l_adj, l_adj_rev);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	std::vector<int> best_scores(n);
	std::vector<edge_t> best_edges(n);
	for (int i : rev_topological_order) {
		for (edge_t& e : l_adj_rev[i]) {
			if (best_scores[e.next] < e.score + best_scores[i]) {
				best_scores[e.next] = e.score + best_scores[i];
				best_edges[e.next] = {i, e.score, e.overlap};
			}
		}
	}

	// 0 is the consensus, we start from there
	edge_t e = best_edges[0];
	std::string ext_consensus = consensus->consensus;
	while (e.overlap) {
		ext_consensus += read_seqs[e.next].substr(e.overlap);
		e = best_edges[e.next];
		consensus->right_ext_reads++;
		if (read_mapqs[e.next] >= config.high_confidence_mapq) consensus->hq_right_ext_reads++;
	}

	if (!consensus->left_clipped) {
		if (consensus->clip_len != consensus_t::UNKNOWN_CLIP_LEN) {
			consensus->clip_len += ext_consensus.length() - consensus->consensus.length();
		}
	} else {
		consensus->end += ext_consensus.length() - consensus->consensus.length();
		if (consensus->max_mapq < config.high_confidence_mapq && consensus->hq_right_ext_reads >= 3) {
			consensus->max_mapq = config.high_confidence_mapq;
		}
	}
	consensus->consensus = ext_consensus;
}

void extend_consensus_to_left(consensus_t* consensus, IntervalTree<ext_read_t*>& candidate_reads_itree,
		hts_pos_t target_start, hts_pos_t target_end, std::string contig_name, hts_pos_t contig_len,
		config_t config, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	std::vector<std::string> read_seqs;
	std::vector<int> read_mapqs;
	read_seqs.push_back(consensus->consensus);
	read_mapqs.push_back(config.high_confidence_mapq);

	get_extension_read_seqs(candidate_reads_itree, read_seqs, read_mapqs, mateseqs_w_mapq, target_start, target_end, contig_len, config, 5000);

	if (read_seqs.size() > 5000) return;

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	build_graph_rev(read_seqs, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, config.read_len/2);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	bool cycle_contains_0 = is_vertex_in_cycle(l_adj, 0);
	if (cycle_contains_0) {
		// remove all edges pointing to 0
		for (std::vector<edge_t>& l_adj_i : l_adj) {
			l_adj_i.erase(std::remove_if(l_adj_i.begin(), l_adj_i.end(), [](edge_t& e) {return e.next == 0;}), l_adj_i.end());
		}
		l_adj_rev[0].clear();
	}

	bool cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) {
		break_cycles(out_edges, l_adj, l_adj_rev);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	std::vector<int> best_scores(n);
	std::vector<edge_t> best_edges(n);
	for (int i : rev_topological_order) {
		for (edge_t& e : l_adj_rev[i]) {
			if (best_scores[e.next] < e.score + best_scores[i]) {
				best_scores[e.next] = e.score + best_scores[i];
				best_edges[e.next] = {i, e.score, e.overlap};
			}
		}
	}

	// 0 is the consensus, we start from there
	edge_t e = best_edges[0];
	std::string ext_consensus = consensus->consensus;
	while (e.overlap) {
		ext_consensus = read_seqs[e.next].substr(0, read_seqs[e.next].length()-e.overlap) + ext_consensus;
		e = best_edges[e.next];
		consensus->left_ext_reads++;
		if (read_mapqs[e.next] >= config.high_confidence_mapq) consensus->hq_left_ext_reads++;
	}

	if (consensus->left_clipped) {
		if (consensus->clip_len != consensus_t::UNKNOWN_CLIP_LEN) {
			consensus->clip_len += ext_consensus.length() - consensus->consensus.length();
		}
	} else {
		consensus->start -= ext_consensus.length() - consensus->consensus.length();
		if (consensus->max_mapq < config.high_confidence_mapq && consensus->hq_left_ext_reads >= 3)
			consensus->max_mapq = config.high_confidence_mapq;
	}
	consensus->consensus = ext_consensus;
}


#endif /* EXTEND_1SR_CONSENSUS_H_ */
