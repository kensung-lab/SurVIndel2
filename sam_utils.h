#ifndef SURVINDEL2_SAM_UTILS_H
#define SURVINDEL2_SAM_UTILS_H

#include <htslib/sam.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

#include "utils.h"

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}

bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}
bool is_duplicate(bam1_t* r) {
	return r->core.flag & BAM_FDUP;
}

bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid;
}
bool is_samestr(bam1_t* r) {
    return is_samechr(r) && (bam_is_rev(r) == bam_is_mrev(r));
}
bool is_outward(bam1_t* r) {
	return is_samechr(r) && ((!bam_is_rev(r) && r->core.isize < 0) || (bam_is_rev(r) && r->core.isize > 0));
}
bool is_long(bam1_t* r, int max_is) {
	return is_samechr(r) && ((!bam_is_rev(r) && r->core.isize > max_is) || (bam_is_rev(r) && r->core.isize < -max_is));
}

int get_left_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' ? bam_cigar_oplen(cigar[0]): 0;
}
int get_right_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' ? bam_cigar_oplen(cigar[r->core.n_cigar-1]): 0;
}

bool is_left_clipped(bam1_t* r, int min_clip_len) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' && bam_cigar_oplen(cigar[0]) >= min_clip_len;
}
bool is_right_clipped(bam1_t* r, int min_clip_len) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' && bam_cigar_oplen(cigar[r->core.n_cigar-1]) >= min_clip_len;
}

bool is_mate_left_clipped(bam1_t* r, int min_clip_len = 0) {
    const uint8_t* mc_tag = bam_aux_get(r, "MC");
    if (mc_tag == NULL) {
        throw "Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.";
    }
    char* mc_tag_str = bam_aux2Z(mc_tag);
    int i = 0, n = 0;
    while (mc_tag_str[i] >= '0' && mc_tag_str[i] <= '9') {
    	n = n*10 + mc_tag_str[i]-'0';
    	i++;
    }
    return mc_tag_str[i] == 'S' && n >= min_clip_len;
}
bool is_mate_right_clipped(bam1_t* r) {
    const uint8_t* mc_tag = bam_aux_get(r, "MC");
    if (mc_tag == NULL) {
        throw "Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.";
    }
    char* mc_tag_str = bam_aux2Z(mc_tag);
    int i = strlen(mc_tag_str)-1;
    return mc_tag_str[i] == 'S';
}

hts_pos_t get_unclipped_start(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return r->core.pos - get_left_clip_size(r);
}
hts_pos_t get_unclipped_end(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_endpos(r) + get_right_clip_size(r);
}

int get_mate_endpos(const bam1_t* r) {
    uint8_t* mcs = bam_aux_get(r, "MC");
    if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

    char* mc = bam_aux2Z(mcs);
    int i = 0, mclen = strlen(mc);

    int len = 0, pos = r->core.mpos;
    while (i < mclen) {
        if (mc[i] >= '0' && mc[i] <= '9') {
            len = (len*10) + (mc[i]-'0');
        } else {
            if (mc[i] != 'I' && mc[i] != 'S') {
                pos += len;
            }
            len = 0;
        }
        i++;
    }
    return pos-1;
}

int64_t get_mq(bam1_t* r) {
    uint8_t* mq = bam_aux_get(r, "MQ");
    if (mq == NULL) {
    	if ((r->core.flag & BAM_FUNMAP) == 0 && (r->core.flag & BAM_FMUNMAP) == 0) {
			std::cerr << "Warning: read pair " << bam_get_qname(r) << " does not have an MQ tag. Please include it." << std::endl;
			std::cerr << (r->core.flag & BAM_FMUNMAP) << std::endl;
    	}
    	return 0;
    }
    return bam_aux2i(mq);
}

void rc(std::string& read) {
    int len = read.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
		char c = std::toupper(read[i]);
		if (c == 'A') read[i] = 'T';
		else if (c == 'C') read[i] = 'G';
		else if (c == 'G') read[i] = 'C';
		else if (c == 'T') read[i] = 'A';
		else c = 'N';
	}
}
void rc(char* read) {
    int len = strlen(read);
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
    	char c = std::toupper(read[i]);
        if (c == 'A') read[i] = 'T';
        else if (c == 'C') read[i] = 'G';
        else if (c == 'G') read[i] = 'C';
        else if (c == 'T') read[i] = 'A';
        else c = 'N';
    }
}

char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
}

std::string get_sequence(bam1_t* r, bool fastq_seq = false) {
    char seq[100000];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (fastq_seq && bam_is_rev(r)) rc(seq);
    return std::string(seq);
}
char* get_sequence_cstr(bam1_t* r, bool fastq_seq = false) {
	char* seq = new char[r->core.l_qseq+1];
	const uint8_t* bam_seq = bam_get_seq(r);
	for (int i = 0; i < r->core.l_qseq; i++) {
		seq[i] = get_base(bam_seq, i);
	}
	seq[r->core.l_qseq] = '\0';
	if (fastq_seq && bam_is_rev(r)) rc(seq);
	return seq;
}

std::string get_cigar_code(bam1_t* r) {
    const uint32_t* cigar = bam_get_cigar(r);
    std::stringstream ss;
    for (int i = 0; i < r->core.n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}

std::pair<int, int> get_mismatches_prefix(bam1_t* read, int prefix_len) {
    int indels = 0, mismatches = 0;
    int original_prefix_len = prefix_len;

    uint32_t* cigar = bam_get_cigar(read);
    for (int j = 0; j < read->core.n_cigar && prefix_len > 0; j++) {
        uint32_t op = cigar[j];
        char opchr = bam_cigar_opchr(op);
        int oplen = bam_cigar_oplen(op);
        if (opchr == 'I' || opchr == 'D') {
            indels += std::min(oplen, prefix_len);
        }
        if (opchr != 'D') {
            prefix_len -= oplen;
        }
    }

    prefix_len = original_prefix_len;
    if (bam_cigar_opchr(cigar[0]) == 'S') prefix_len -= bam_cigar_oplen(cigar[0]) ;

    char* md_tag = bam_aux2Z(bam_aux_get(read, "MD"));
    int md_tag_len = strlen(md_tag);
    int n = 0;
    for (int i = 0; i < md_tag_len && prefix_len > 0; i++) {
        if (md_tag[i] >= '0' && md_tag[i] <= '9') {
            n = n*10 + md_tag[i]-'0';
        } else {
            prefix_len -= n;
            if (prefix_len <= 0) break;
            n = 0;

            if (md_tag[i] == '^') {
                while (md_tag[i] < '0' || md_tag[i] > '9') i++;
            } else {
                mismatches++;
                prefix_len--;
            }
        }
    }

    return {mismatches, indels};
}

std::pair<int, int> get_mismatches_suffix(bam1_t* read, int prefix_len) {
    int indels = 0, mismatches = 0;
    int original_prefix_len = prefix_len;

    uint32_t* cigar = bam_get_cigar(read);
    for (int j = read->core.n_cigar-1; j >= 0 && prefix_len > 0; j--) {
        uint32_t op = cigar[j];
        char opchr = bam_cigar_opchr(op);
        int oplen = bam_cigar_oplen(op);
        if (opchr == 'I' || opchr == 'D') {
            indels += std::min(oplen, prefix_len);
        }
        if (opchr != 'D') {
            prefix_len -= oplen;
        }
    }

    prefix_len = original_prefix_len;
    if (bam_cigar_opchr(cigar[0]) == 'S') prefix_len -= bam_cigar_oplen(cigar[0]) ;

    char* md_tag = bam_aux2Z(bam_aux_get(read, "MD"));
    int md_tag_len = strlen(md_tag);
    int n = 0, m = 1;
    for (int i = md_tag_len-1; i >= 0 && prefix_len > 0; i--) {
        if (md_tag[i] >= '0' && md_tag[i] <= '9') {
            n += m * (md_tag[i]-'0');
            m *= 10;
        } else {
            prefix_len -= n;
            if (prefix_len <= 0) break;
            n = 0, m = 1;

            if (i == 0 || (md_tag[i] != '^' && md_tag[i-1] >= '0' && md_tag[i-1] <= '9')) {
                mismatches++;
            }
        }
    }

    return {mismatches, indels};
}

int64_t get_AS_tag(bam1_t* read) {
    uint8_t* aux_get = bam_aux_get(read, "AS");
    return bam_aux2i(aux_get);
}
int get_aligned_portion_len(bam1_t* read) {
    return read->core.l_qseq - get_left_clip_size(read) - get_right_clip_size(read);
}

samFile* get_writer(std::string path, bam_hdr_t* header) {
    samFile* writer = sam_open(path.c_str(), "wb");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + path;
    }
    return writer;
}

struct open_samFile_t {
    samFile* file;
    bam_hdr_t* header;
    hts_idx_t* idx;

    open_samFile_t() {}
};

open_samFile_t* open_samFile(std::string fname_str, bool index_file = false) {
    const char* fname = fname_str.c_str();
    open_samFile_t* sam_file = new open_samFile_t;
    sam_file->file = sam_open(fname, "r");
    if (sam_file->file == NULL) {
        throw "Could not open " + std::string(fname);
    }

    if (index_file) {
        int code = sam_index_build(fname, 0);
        if (code != 0) {
            throw "Cannot index " + std::string(fname);
        }
    }

    sam_file->idx = sam_index_load(sam_file->file, sam_file->file->fn);
    if (sam_file->idx == NULL) {
        throw "Unable to open index for " + std::string(fname);
    }

    sam_file->header = sam_hdr_read(sam_file->file);
    if (sam_file->header == NULL) {
        throw "Unable to open header for " + std::string(fname);
    }

    return sam_file;
}

void close_samFile(open_samFile_t* f) {
    hts_idx_destroy(f->idx);
    bam_hdr_destroy(f->header);
    sam_close(f->file);
    delete f;
}

struct bam_redux_t {
    static const uint8_t IS_REV = 1, IS_MREV = 2, IS_INTER_CHR = 4, IS_MATE_LC = 8, IS_MATE_RC = 16;

    hts_pos_t start, end, mstart, isize;
    int left_clip_size, right_clip_size;
    int nm = 0;
    uint8_t flag = 0, mapq = 0;
    std::vector<uint8_t> seq;
    std::vector<uint8_t> qual;
    std::vector<uint32_t> cigar;

    // TODO: temporary
    std::string qname;
    int as = 0;

    bam_redux_t() {}
    bam_redux_t(bam1_t* read) : start(read->core.pos), end(bam_endpos(read)), mstart(read->core.mpos),
        isize(read->core.isize), mapq(read->core.qual), left_clip_size(get_left_clip_size(read)), right_clip_size(get_right_clip_size(read)),
		nm(bam_aux2i(bam_aux_get(read, "NM"))), as(bam_aux2i(bam_aux_get(read, "AS"))) {

    	qname = bam_get_qname(read);
        if (bam_is_rev(read)) flag |= IS_REV;
        if (bam_is_mrev(read)) flag |= IS_MREV;
        if (read->core.tid != read->core.mtid) flag |= IS_INTER_CHR;
        if (is_mate_left_clipped(read)) flag |= IS_MATE_LC;
        if (is_mate_right_clipped(read)) flag |= IS_MATE_RC;

        uint8_t* seq_array = bam_get_seq(read);
        seq = std::vector<uint8_t>(seq_array, seq_array+(read->core.l_qseq+1)/2);

        uint8_t* qual_array = bam_get_qual(read);
        qual = std::vector<uint8_t>(qual_array, qual_array+read->core.l_qseq);

        uint32_t* cigar_array = bam_get_cigar(read);
        cigar = std::vector<uint32_t>(cigar_array, cigar_array+read->core.n_cigar);
    }

    int seq_len() {
        return qual.size();
    }

    bool is_rev() {
        return flag & IS_REV;
    }
    bool is_mrev() {
        return flag & IS_MREV;
    }
    bool is_inter_chr() {
        return flag & IS_INTER_CHR;
    }
    bool mate_left_clipped() {
        return flag & IS_MATE_LC;
    }
    bool mate_right_clipped() {
        return flag & IS_MATE_RC;
    }

    hts_pos_t unclipped_start() {
        return start-left_clip_size;
    }
    hts_pos_t unclipped_end() {
        return end+right_clip_size;
    }

    std::string cigar_string() {
        std::stringstream ss;
        for (uint32_t c : cigar) ss << bam_cigar_oplen(c) << bam_cigar_opchr(c);
        return ss.str();
    }
};

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

#endif //SURVINDEL2_SAM_UTILS_H
