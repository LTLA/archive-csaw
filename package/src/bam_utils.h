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
    
struct BamRead {
    BamRead();
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

void decompose_cigar(bam1_t*, int&, int&);

#endif
