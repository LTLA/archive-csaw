#include <string>
#include <map>

#include "csaw.h"
#include "sam.h"
#include "bgzf.h"

extern "C" {

class Bamfile {
public:
    Bamfile(const char * path, const char* xpath) {
        in = sam_open(path, "rb");
        if (in == NULL) {
            throw std::runtime_error("failed to open BAM file path");
        }
        try {
            index = bam_index_load(xpath); 
            if (index==NULL) { throw std::runtime_error("failed to open BAM index file path"); }
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
        read=bam_init1();
        return;
    }

    ~Bamfile() {
        bam_hdr_destroy(header);
        hts_idx_destroy(index); 
        sam_close(in);
        bam_destroy1(read);
    }

    samFile* in;
    bam1_t* read;
    hts_idx_t* index;
    bam_hdr_t * header;
};

class BamIterator {
public:
    BamIterator(const Bamfile& bf, const char* chr) : iter(NULL) {
        int cid=bam_name2id(bf.header, chr);
        iter=bam_itr_queryi(bf.index, cid, 0, (bf.header->target_len)[cid]); 
    }
    ~BamIterator() { 
        bam_itr_destroy(iter); 
    }
    hts_itr_t* iter;
};

bool is_mapped(bam1_t* read, bool use_qual, int minqual, bool rmdup) {
    if (use_qual && (read -> core).qual < minqual) { return false; }
    if (rmdup && ((read -> core).flag & BAM_FDUP)==0) { return false; }
    return true;
}

/* Strolls through the file for each chromosome and accumulates paired-end statistics; 
 * forward and reverse reads (position and width), singles, unoriented, and names of
 * inter-chromosomals.
 */

SEXP extract_pair_data(SEXP bam, SEXP index, SEXP chr, SEXP mapq, SEXP dedup, SEXP get_names) try {
    // Checking input values.
    if (!isString(bam) || LENGTH(bam)!=1) {
        throw std::runtime_error("BAM file path must be a string");
    }
    if (!isString(index) || LENGTH(index)!=1) {
        throw std::runtime_error("BAM index file path must be a string"); 
    }
    if (!isString(chr) || LENGTH(chr)!=1) { 
        throw std::runtime_error("chromosome name should be a string"); 
    }

    if (!isInteger(mapq) || LENGTH(mapq)!=1) {
        throw std::runtime_error("mapping quality should be an integer scalar");
    }    
    const int minqual=asInteger(mapq);
    const bool use_qual=(minqual==NA_INTEGER);

    if (!isLogical(dedup) || LENGTH(dedup)!=1) {
        throw std::runtime_error("duplicate removal should be a logical scalar"); 
    }
    const bool rmdup=asLogical(dedup);

    if (!isLogical(get_names) || LENGTH(get_names)!=1) { 
        throw std::runtime_error("get-stat specification should be a logical scalar"); 
    }
    const bool getnames=asLogical(get_names);

    // Initializing odds and ends.
    Bamfile bf(CHAR(STRING_ELT(bam, 0)), CHAR(STRING_ELT(index, 0)));
    BamIterator biter(bf, CHAR(STRING_ELT(chr, 0)));
        
    typedef std::map<std::pair<int, std::string>, std::pair<int, int> > Holder;
    Holder first_holder, second_holder;
    std::deque<std::string> interchr_names_1, interchr_names_2;
    int singles=0, both_mapped=0, one_unmapped=0, totals=0, unoriented=0;

    std::pair<int, std::string> current;
    Holder::iterator ith;
    std::deque<int> forward_pos_out, forward_len_out, reverse_pos_out, reverse_len_out;
    int forward_pos, forward_len, reverse_pos, reverse_len, curpos, curlen;
    bool am_mapped, is_first;

    while (bam_itr_next(bf.in, biter.iter, bf.read) >= 0){
        ++totals;

        /* Reasons to not add a read: */
       
        // If it's a singleton.
        if (((bf.read -> core).flag & BAM_FPAIRED)==0) {
            ++singles;
            continue;
        }

        // Or if it's inter-chromosomal.
        am_mapped=is_mapped(bf.read, use_qual, minqual, rmdup);
        is_first=(((bf.read->core).flag & BAM_FREAD1)!=0);
        if (is_first==(((bf.read->core).flag & BAM_FREAD2)!=0)) { 
            throw std::runtime_error("exactly one of the first/second BAM fields must be set"); 
        }
        
        if ((bf.read -> core).mtid!=(bf.read -> core).tid) { 
            if (getnames && am_mapped) { 
                (is_first ? interchr_names_1 : interchr_names_2).push_back(std::string(bam_get_qname(bf.read))); 
            } 
            continue;
        }

        // Or, if we can see that it is unmapped (diagnostics depend on partner).
        if (((bf.read->core).flag & BAM_FUNMAP) >= 0) { 
            continue;
        } 
        
        // Or, if we can see that its partner is obviously unmapped.
        if (((bf.read -> core).flag & BAM_FMUNMAP) >= 0) {
            if (am_mapped) { ++one_unmapped; }
            continue;
        }

        /* Checking the map and adding it if it doesn't exist. */

        curpos = (bf.read -> core).pos;
        curlen = (bf.read -> core).l_qseq;

        if ((bf.read -> core).mpos < (bf.read -> core).pos) {
            current.first = (bf.read -> core).mpos;
            current.second.assign(bam_get_qname(bf.read));
            Holder& holder=( is_first ? second_holder : first_holder );
            ith=holder.find(current);

            if (ith != holder.end()) {
                if (!am_mapped) {
                    ++one_unmapped;
                    holder.erase(ith);
                    continue;
                }
                ++both_mapped;
            
                if ( ((ith->second).second < 0) == bool(bam_is_rev(bf.read)) ) {
                    ++unoriented; 
                    holder.erase(ith);
                    continue;
                } 

                if ( ((ith->second).second < 0) ) { 
                    forward_pos = curpos;
                    forward_len = curlen;
                    reverse_pos = (ith->second).first;
                    reverse_len = -((ith->second).second);
                } else {
                    forward_pos = (ith->second).first;
                    forward_len = ((ith->second).second);
                    reverse_pos = curpos;
                    reverse_len = curlen;
                }

                if (forward_pos > reverse_pos || forward_pos+forward_len > reverse_pos + reverse_len) {
                    ++unoriented;
                    continue;
                }

                forward_pos_out.push_back(forward_pos);
                forward_len_out.push_back(forward_pos);
                reverse_pos_out.push_back(reverse_pos);
                reverse_len_out.push_back(reverse_pos);
                holder.erase(ith);
            } else if (am_mapped) {
                ++one_unmapped;
            }
        } else if (am_mapped) {
            current.first = (bf.read -> core).pos;
            current.second.assign(bam_get_qname(bf.read));
            Holder& holder=( is_first ? first_holder : second_holder );
            holder[current] = std::make_pair(curpos, curlen * (bam_is_rev(bf.read) ? -1 : 1));
        }
    }

    // Leftovers treated as one_unmapped.
    one_unmapped += first_holder.size() + second_holder.size();
    first_holder.clear();
    second_holder.clear();

    SEXP output=PROTECT(allocVector(VECSXP, 4));
    try {
        for (size_t left=0; left<1; ++left) { 
            SET_VECTOR_ELT(output, left, allocVector(VECSXP, 2));
            SEXP curout=VECTOR_ELT(output, left);
            std::deque<int>& curposout=(left==0 ? forward_pos_out : reverse_pos_out);
            SET_VECTOR_ELT(curout, 0, allocVector(INTSXP, curposout.size()));
            std::copy(curposout.begin(), curposout.end(), INTEGER(VECTOR_ELT(curout, 0)));
            std::deque<int>& curlenout=(left==0 ? forward_len_out : reverse_len_out);
            SET_VECTOR_ELT(curout, 1, allocVector(INTSXP, curlenout.size()));
            std::copy(curlenout.begin(), curlenout.end(), INTEGER(VECTOR_ELT(curout, 1)));
        }

        SET_VECTOR_ELT(output, 2, allocVector(INTSXP, 5));
        int * diagptr=INTEGER(VECTOR_ELT(output, 2));
        diagptr[0] = totals;
        diagptr[1] = both_mapped;
        diagptr[2] = singles;
        diagptr[3] = one_unmapped;
        diagptr[4] = unoriented;

        SET_VECTOR_ELT(output, 3, allocVector(VECSXP, 2));
        SEXP allnames=VECTOR_ELT(output, 3);
        SET_VECTOR_ELT(allnames, 0, allocVector(STRSXP, interchr_names_1.size()));
        SET_VECTOR_ELT(allnames, 1, allocVector(STRSXP, interchr_names_2.size()));
        SEXP outnames=VECTOR_ELT(allnames, 0);
        for (size_t i=0; i<interchr_names_1.size(); ++i) {
            SET_STRING_ELT(outnames, i, mkChar(interchr_names_1[i].c_str()));
        }
        outnames=VECTOR_ELT(allnames, 1);
        for (size_t i=0; i<interchr_names_2.size(); ++i) {
            SET_STRING_ELT(outnames, i, mkChar(interchr_names_2[i].c_str()));
        }
    } catch (std::exception &e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception &e) {
    return mkString(e.what());
}


}
