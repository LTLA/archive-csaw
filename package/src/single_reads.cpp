#include "bam_utils.h"

SEXP extract_single_data(SEXP bam, SEXP index, SEXP chr, SEXP start, SEXP end, 
        SEXP mapq, SEXP dedup, SEXP use_forward, SEXP use_first) try {
    // Checking input values.
    if (!isInteger(mapq) || LENGTH(mapq)!=1) {
        throw std::runtime_error("mapping quality should be an integer scalar");
    }    
    const int minqual=asInteger(mapq);

    if (!isLogical(dedup) || LENGTH(dedup)!=1) {
        throw std::runtime_error("duplicate removal should be a logical scalar"); 
    }
    const bool rmdup=asLogical(dedup);

    int set_flags=0, unset_flags=0;
    if (!isLogical(use_forward) || LENGTH(use_forward)!=1) {
        throw std::runtime_error("forward usage specification should be a logical scalar"); 
    } else {
        const int useforward=asLogical(use_forward);
        if (useforward!=NA_LOGICAL) {
            if (useforward) {
                unset_flags+=BAM_FREVERSE;                
            } else {
                set_flags+=BAM_FREVERSE;                
            }
        }
    }

    if (!isLogical(use_first) || LENGTH(use_first)!=1) {
        throw std::runtime_error("first usage specification should be a logical scalar"); 
    } else {
        const int usefirst=asLogical(use_first);
        if (usefirst!=NA_LOGICAL) {
            if (usefirst) {
                set_flags+=BAM_FREAD1;
                unset_flags+=BAM_FREAD2;
            } else {
                set_flags+=BAM_FREAD2;
                unset_flags+=BAM_FREAD1;
            }
        }
    }

    // Initializing odds and ends.
    BamFile bf(bam, index);
    BamRead br;
    BamIterator biter(bf, chr, start, end);
    std::deque<int> forward_pos, forward_len, forward_off,
        reverse_pos, reverse_len, reverse_off;
    int curpos, curlen, curoff;

    while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){    
//        // If we can see that it is obviously unmapped (IMPOSSIBLE for a sorted file).
//        if (((br.read -> core).flag & BAM_FUNMAP)!=0) { 
//            continue;
//        } 
        
        if (!br.is_well_mapped(minqual, rmdup)) { continue; }
        if (((br.read->core).flag & set_flags)!=set_flags) { continue; }
        if (((br.read->core).flag & unset_flags)!=0) { continue; }
        curpos = (br.read -> core).pos + 1;
        br.decompose_cigar(curlen, curoff);
        
        if (bam_is_rev(br.read)) { 
            reverse_pos.push_back(curpos);
            reverse_len.push_back(curlen);
            reverse_off.push_back(curoff);
        } else {
            forward_pos.push_back(curpos);
            forward_len.push_back(curlen);
            forward_off.push_back(curoff);
        }        
    }

    // Storing all output.
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, 3));
        SEXP left=VECTOR_ELT(output, 0);
        store_int_output(left, 0, forward_pos);
        store_int_output(left, 1, forward_len);
        store_int_output(left, 2, forward_off);
        
        SET_VECTOR_ELT(output, 1, allocVector(VECSXP, 3));
        SEXP right=VECTOR_ELT(output, 1);
        store_int_output(right, 0, reverse_pos);
        store_int_output(right, 1, reverse_len);
        store_int_output(right, 2, reverse_off);
    
    } catch (std::exception &e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception &e) {
    return mkString(e.what());
}

