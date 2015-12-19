#include "bam_utils.h"

BamFile::BamFile(SEXP bam, SEXP idx) {
    if (!isString(bam) || LENGTH(bam)!=1) {
        throw std::runtime_error("BAM file path must be a string");
    }
    const char* path=CHAR(STRING_ELT(bam, 0));
    if (!isString(idx) || LENGTH(idx)!=1) {
        throw std::runtime_error("BAM index file path must be a string"); 
    }
    const char* xpath=CHAR(STRING_ELT(idx, 0));

    in = sam_open(path, "rb");
    if (in == NULL) {
        std::stringstream err;
        err << "failed to open BAM file at '" << path << "'";
        throw std::runtime_error(err.str());
    }
    try {
        index = bam_index_load(xpath); 
        if (index==NULL) { 
            std::stringstream err;
            err << "failed to open BAM index at '" << xpath << "'";
            throw std::runtime_error(err.str());
        }
        try {
            header=sam_hdr_read(in);
        } catch (std::exception& e) {
            hts_idx_destroy(index);
            throw;
        }
    } catch (std::exception& e) {
        sam_close(in);
        throw;
    }
    bgzf_set_cache_size(in->fp.bgzf, 100*BGZF_MAX_BLOCK_SIZE);
    return;
}

BamFile::~BamFile() {
    bam_hdr_destroy(header);
    hts_idx_destroy(index); 
    sam_close(in);
    return;
}

BamRead::BamRead() {
    read=bam_init1();
    return;
}

bool BamRead::is_well_mapped(const int& minqual, const bool& rmdup) const {
    if (minqual!=NA_INTEGER && (read -> core).qual < minqual) { return false; } 
    if (rmdup && ((read -> core).flag & BAM_FDUP)!=0) { return false; }
    return true;
}

void BamRead::decompose_cigar(int& alen, int& offset) const {
    const int n_cigar=(read->core).n_cigar;
    if (n_cigar==0) { 
        std::stringstream err;
        err << "zero-length CIGAR for mapped read '" << bam_get_qname(read) << "'";
        throw std::runtime_error(err.str());
    }
    uint32_t* cigar=bam_get_cigar(read);
    
    alen=bam_cigar2rlen(n_cigar, cigar);
    offset=0;
    if (bam_is_rev(read)) {
        if (bam_cigar_op(cigar[n_cigar-1])==BAM_CSOFT_CLIP) { offset = bam_cigar_oplen(cigar[n_cigar-1]); }
    } else {
        if (bam_cigar_op(cigar[0])==BAM_CSOFT_CLIP) { offset = bam_cigar_oplen(cigar[0]); }
    }
    return;
}

BamRead::~BamRead() { 
    bam_destroy1(read);
    return;
}

BamIterator::BamIterator(const BamFile& bf) : iter(NULL) {
    iter=bam_itr_queryi(bf.index, HTS_IDX_NOCOOR, 0, 0);
    return;
}

BamIterator::BamIterator(const BamFile& bf, SEXP Chr, SEXP Start, SEXP End) : iter(NULL) {
    if (!isString(Chr) || LENGTH(Chr)!=1) { 
        throw std::runtime_error("chromosome name should be a string"); 
    }
    const char* chr=CHAR(STRING_ELT(Chr, 0));
    if (!isInteger(Start) || LENGTH(Start)!=1) { 
        throw std::runtime_error("region start should be an integer scalar"); 
    }
    int start=asInteger(Start)-1;
    if (!isInteger(End) || LENGTH(End)!=1) { 
        throw std::runtime_error("region end should be an integer scalar"); 
    }
    int end=asInteger(End);

    int cid=bam_name2id(bf.header, chr);
    if (cid==-1) {
        std::stringstream err;
        err << "reference sequence '" << chr << "' missing in BAM header";
        throw std::runtime_error(err.str());
    }

    if (start < 0) { start=0; }
    const int curlen = (bf.header->target_len)[cid];
    if (end > curlen) { end = curlen; }
    if (start > end) {
        throw std::runtime_error("invalid values for region start/end coordinates");
    }
    iter=bam_itr_queryi(bf.index, cid, start, end);
}

BamIterator::BamIterator(const BamFile& bf, int cid, int start, int end) : iter(NULL) {
    iter=bam_itr_queryi(bf.index, cid, start, end);
}

BamIterator::~BamIterator() { 
    bam_itr_destroy(iter); 
}

void store_int_output(SEXP& dest, int index, const std::deque<int>& host) {
    SET_VECTOR_ELT(dest, index, allocVector(INTSXP, host.size()));
    std::copy(host.begin(), host.end(), INTEGER(VECTOR_ELT(dest, index)));
    return;
}

void store_names(SEXP& dest, int index, const std::deque<std::string>& names) {
    SET_VECTOR_ELT(dest, index, allocVector(STRSXP, names.size()));
    SEXP current=VECTOR_ELT(dest, index);
    for (size_t i=0; i<names.size(); ++i) {
        SET_STRING_ELT(current, i, mkChar(names[i].c_str()));
    }
    return;
}

