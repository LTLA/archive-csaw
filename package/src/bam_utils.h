#ifndef BAM_UTILS_H
#define BAM_UTILS_H

#include "csaw.h"
#include "sam.h"
#include "hts.h"
#include "bgzf.h"

struct BamFile {
    BamFile(SEXP, SEXP);
    ~BamFile();
    samFile* in;
    hts_idx_t* index;
    bam_hdr_t * header;
};

struct AlignData {
    AlignData();
    int len;
    bool is_reverse;
};
    
struct BamRead {
    BamRead();
    bool is_well_mapped(const int&, const bool&) const;
    void extract_data(AlignData&) const;
    ~BamRead();
    bam1_t* read;
};

struct BamIterator {
    BamIterator(const BamFile&);
    BamIterator(const BamFile&, SEXP, SEXP, SEXP);
    BamIterator(const BamFile&, int, int, int);
    ~BamIterator();
    hts_itr_t* iter;
};

void store_int_output(SEXP&, int, const std::deque<int>&);

void store_names(SEXP&, int, const std::deque<std::string>&);

#endif
