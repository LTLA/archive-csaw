#include <string>
#include <sstream>
#include <map>

#include "csaw.h"
#include "sam.h"
#include "bgzf.h"

extern "C" {

class Bamfile {
public:
    Bamfile(SEXP bam, SEXP idx) {
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
    BamIterator(const Bamfile& bf, SEXP Chr, SEXP Start, SEXP End) : iter(NULL) {
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

void decompose_cigar(bam1_t* read, int& alen, int& offset) {
    const int n_cigar=(read->core).n_cigar;
    if (n_cigar==0) { 
        std::stringstream err;
        err << "zero-length CIGAR for read '" << bam_get_qname(read) << "'";
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

class OutputContainer {
public:
    OutputContainer(bool d) : diagnostics(d), totals(0) {}

    void add_genuine(int pos1, int len1, int off1, int pos2, int len2, int off2, bool isreverse1, bool isfirst1) {
        mate_reverse = (len2 < 0);
        if (mate_reverse) { len2*=-1; }
        if (mate_reverse==isreverse1) {
            add_unoriented(pos1, len1, pos2, len2, isfirst1); 
            return;
        }
         
        if (mate_reverse) { 
            forward_pos = pos1;
            forward_len = len1;
            forward_off = off1;
            reverse_pos = pos2;
            reverse_len = len2;
            reverse_off = off2;
        } else {
            forward_pos = pos2;
            forward_len = len2;
            forward_off = off2;
            reverse_pos = pos1;
            reverse_len = len1;
            reverse_off = off1;
        }

        if (forward_pos > reverse_pos || forward_pos+forward_len > reverse_pos + reverse_len) {
            add_unoriented(pos1, len1, pos2, len2, isfirst1);
            return; 
        }

        ++forward_pos; // Get to 1-indexed positions.
        ++reverse_pos;
        forward_pos_out.push_back(forward_pos);
        forward_len_out.push_back(forward_len);
        forward_off_out.push_back(forward_off);
        reverse_pos_out.push_back(reverse_pos);
        reverse_len_out.push_back(reverse_len);
        reverse_off_out.push_back(reverse_off);
        return;
    }

    void add_unoriented(int pos1, int len1, int pos2, int len2, bool isfirst1) {
        if (!diagnostics) { return; }
        ++pos1;
        ++pos2;
        if (isfirst1) {
            ufirst_pos.push_back(pos1);
            ufirst_len.push_back(len1);
            usecond_pos.push_back(pos2);
            usecond_len.push_back(len2);
        } else {
            ufirst_pos.push_back(pos2);
            ufirst_len.push_back(len2);
            usecond_pos.push_back(pos1);
            usecond_len.push_back(len1);
        }
    }

    void add_onemapped(int pos, int len) {
        if (!diagnostics) { return; }
        ++pos;
        if (len < 0) { len*=-1; }
        onemap_pos.push_back(pos);
        onemap_len.push_back(len);
        return;
    }

    void add_single(int pos, int len) {
        if (!diagnostics) { return; }
        ++pos;
        single_pos.push_back(pos);
        single_len.push_back(len);
        return;
    }

    void add_interchr(int pos, int len, const char* name, bool isfirst) {
        if (!diagnostics) { return; }
        ++pos;
        if (isfirst) { 
            ifirst_pos.push_back(pos);
            ifirst_len.push_back(len);
            interchr_names_1.push_back(std::string(name));
        } else {
            isecond_pos.push_back(pos);
            isecond_len.push_back(len);
            interchr_names_2.push_back(std::string(name));
        }
        return;
    }

    void store_output(SEXP& dest, int index, const std::deque<int>& host) {
        SET_VECTOR_ELT(dest, index, allocVector(INTSXP, host.size()));
        std::copy(host.begin(), host.end(), INTEGER(VECTOR_ELT(dest, index)));
        return;
    }

    void store_names(SEXP dest, int index, std::deque<std::string>& names) {
        SET_VECTOR_ELT(dest, index, allocVector(STRSXP, names.size()));
        SEXP current=VECTOR_ELT(dest, index);
        for (size_t i=0; i<names.size(); ++i) {
            SET_STRING_ELT(current, i, mkChar(names[i].c_str()));
        }
        return;
    }

    const bool diagnostics;
    int totals;
    bool mate_reverse;
    int forward_pos, forward_len, forward_off, reverse_pos, reverse_len, reverse_off;
    std::deque<int> forward_pos_out, forward_len_out, forward_off_out, 
        reverse_pos_out, reverse_len_out, reverse_off_out;
    std::deque<int> ufirst_pos, ufirst_len, usecond_pos, usecond_len;
    std::deque<int> onemap_pos, onemap_len;
    std::deque<int> single_pos, single_len;
    std::deque<std::string> interchr_names_1, interchr_names_2;
    std::deque<int> ifirst_pos, ifirst_len, isecond_pos, isecond_len;
};


/* Strolls through the file for each chromosome and accumulates paired-end statistics; 
 * forward and reverse reads (position and width), singles, unoriented, and names of
 * inter-chromosomals.
 */

SEXP extract_pair_data(SEXP bam, SEXP index, SEXP chr, SEXP start, SEXP end, SEXP mapq, SEXP dedup, SEXP get_names) try {
    // Checking input values.
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
    Bamfile bf(bam, index);
    BamIterator biter(bf, chr, start, end);
    OutputContainer oc(getnames);
        
    typedef std::map<std::pair<int, std::string>, std::pair<int, int> > Holder;
    Holder first_holder, second_holder;
    std::pair<int, std::string> current;
    Holder::iterator ith;
    int forward_pos, forward_len, reverse_pos, reverse_len, curpos, curlen, curoff;
    bool am_mapped, is_first;

    bool mate_is_in;
    std::set<std::string> identical_pos;
    std::set<std::string>::iterator itip;
    int last_identipos=-1;

    while (bam_itr_next(bf.in, biter.iter, bf.read) >= 0){
        ++oc.totals;

        /* Reasons to not add a read: */
       
        // If we can see that it is obviously unmapped.
        if (((bf.read -> core).flag & BAM_FUNMAP)!=0) { 
            // We don't filter by additional mapping criteria, as we need to search 'holder' to pop out the partner and to store diagnostics.
            continue;
        } 
        
        // (just getting some stats here).
        curpos = (bf.read -> core).pos;
        decompose_cigar(bf.read, curlen, curoff);
        am_mapped=is_mapped(bf.read, use_qual, minqual, rmdup);

        // Or If it's a singleton.
        if (((bf.read -> core).flag & BAM_FPAIRED)==0) {
            if (am_mapped) { oc.add_single(curpos, curlen); }
            continue;
        }

        // Or, if we can see that its partner is obviously unmapped.
        if (((bf.read -> core).flag & BAM_FMUNMAP)!=0) {
            if (am_mapped) { oc.add_onemapped(curpos, curlen); }
            continue;
        }

        // Or if it's inter-chromosomal.
        is_first=(((bf.read->core).flag & BAM_FREAD1)!=0);
        if (is_first==(((bf.read->core).flag & BAM_FREAD2)!=0)) { 
            std::stringstream err;
            err << "read '" << bam_get_qname(bf.read) << "' must be either first or second in the pair";
            throw std::runtime_error(err.str()); 
        }
      
        if ((bf.read -> core).mtid!=(bf.read -> core).tid) { 
            if (am_mapped) { oc.add_interchr(curpos, curlen, bam_get_qname(bf.read), is_first); } 
            continue;
        }

        /* Checking the map and adding it if it doesn't exist. */
        
        current.second.assign(bam_get_qname(bf.read));
        mate_is_in=false;
        if ((bf.read -> core).mpos < curpos) {
            mate_is_in=true;
        } else if ((bf.read -> core).mpos == curpos) {
            // Identical mpos to curpos needs careful handling to figure out whether we've already seen it.
            if (curpos!=last_identipos) { 
                identical_pos.clear();
                last_identipos=curpos;
            }
            itip=identical_pos.lower_bound(current.second);
            if (itip!=identical_pos.end() && !(identical_pos.key_comp()(current.second, *itip))) {
                mate_is_in=true;
                identical_pos.erase(itip);
            } else {
                identical_pos.insert(itip, current.second);
            }
        }

        if (mate_is_in) {
            current.first = (bf.read -> core).mpos;
            Holder& holder=( is_first ? second_holder : first_holder );
            ith=holder.find(current);

            if (ith != holder.end()) { 
                if (!am_mapped) {
                    // Searching to pop out the mate, to reduce the size of 'holder' for the remaining searches (and to store diagnostics).
                    oc.add_onemapped((ith->first).first, (ith->second).first);
                    holder.erase(ith);
                    continue;
                }

                oc.add_genuine(curpos, curlen, curoff, 
                        (ith->first).first, (ith->second).first, (ith->second).second, 
                        bool(bam_is_rev(bf.read)), is_first);
                holder.erase(ith);
            } else if (am_mapped) {
                // Only possible if the mate didn't get added because 'am_mapped' was false.
                oc.add_onemapped(curpos, curlen);
            }
        } else if (am_mapped) {
            current.first = curpos;
            Holder& holder=( is_first ? first_holder : second_holder );
            holder[current] = std::make_pair(curlen * (bam_is_rev(bf.read) ? -1 : 1), curoff);
        }
    }

    // Leftovers treated as one_unmapped; marked as paired, but the mate is not in file.
    for (ith=first_holder.begin(); ith!=first_holder.end(); ++ith) { 
        oc.add_onemapped((ith->first).first, (ith->second).first);
    }
    for (ith=second_holder.begin(); ith!=second_holder.end(); ++ith) { 
        oc.add_onemapped((ith->first).first, (ith->second).first);
    }  
    first_holder.clear();
    second_holder.clear();

    // Storing all output.
    SEXP output=PROTECT(allocVector(VECSXP, getnames ? 9 : 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, 3));
        SEXP left=VECTOR_ELT(output, 0);
        oc.store_output(left, 0, oc.forward_pos_out);
        oc.store_output(left, 1, oc.forward_len_out);
        oc.store_output(left, 2, oc.forward_off_out);
        
        SET_VECTOR_ELT(output, 1, allocVector(VECSXP, 3));
        SEXP right=VECTOR_ELT(output, 1);
        oc.store_output(right, 0, oc.reverse_pos_out);
        oc.store_output(right, 1, oc.reverse_len_out);
        oc.store_output(right, 2, oc.reverse_off_out);
    
        if (getnames) {
            SET_VECTOR_ELT(output, 2, ScalarInteger(oc.totals));
            
            SET_VECTOR_ELT(output, 3, allocVector(VECSXP, 2));
            SEXP singles=VECTOR_ELT(output, 3);
            oc.store_output(singles, 0, oc.single_pos);
            oc.store_output(singles, 1, oc.single_len);

            SET_VECTOR_ELT(output, 4, allocVector(VECSXP, 2));
            SEXP first=VECTOR_ELT(output, 4);
            oc.store_output(first, 0, oc.ufirst_pos);
            oc.store_output(first, 1, oc.ufirst_len);
            
            SET_VECTOR_ELT(output, 5, allocVector(VECSXP, 2));
            SEXP second=VECTOR_ELT(output, 5);
            oc.store_output(second, 0, oc.usecond_pos);
            oc.store_output(second, 1, oc.usecond_len);

            SET_VECTOR_ELT(output, 6, allocVector(VECSXP, 2));
            SEXP onemap=VECTOR_ELT(output, 6);
            oc.store_output(onemap, 0, oc.onemap_pos);
            oc.store_output(onemap, 1, oc.onemap_len);

            SET_VECTOR_ELT(output, 7, allocVector(VECSXP, 3));
            SEXP interchr1=VECTOR_ELT(output, 7);
            oc.store_output(interchr1, 0, oc.ifirst_pos);
            oc.store_output(interchr1, 1, oc.ifirst_len);
            oc.store_names(interchr1, 2, oc.interchr_names_1);

            SET_VECTOR_ELT(output, 8, allocVector(VECSXP, 3));
            SEXP interchr2=VECTOR_ELT(output, 8);
            oc.store_output(interchr2, 0, oc.isecond_pos);
            oc.store_output(interchr2, 1, oc.isecond_len);
            oc.store_names(interchr2, 2, oc.interchr_names_2);
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
