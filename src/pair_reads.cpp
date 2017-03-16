#include "bam_utils.h"

struct OutputContainer {
    OutputContainer(bool d) : diagnostics(d), totals(0) {}

    void add_genuine(int pos1, const AlignData& data1, int pos2, const AlignData& data2, bool isfirst1) {
        if (data2.is_reverse==data1.is_reverse) {
            add_unoriented(pos1, data1, pos2, data2, isfirst1); 
            return;
        }
         
        if (data2.is_reverse) {
            forward_pos = pos1;
            forward_len = data1.len;
            reverse_pos = pos2;
            reverse_len = data2.len;
        } else {
            forward_pos = pos2;
            forward_len = data2.len;
            reverse_pos = pos1;
            reverse_len = data1.len;
        }

        // Completely overrun fragments go into the 'unoriented' basket.
        if (forward_pos >= reverse_pos + reverse_len) {
            add_unoriented(pos1, data1, pos2, data2, isfirst1);
            return; 
        } 

        forward_pos_out.push_back(forward_pos); 
        forward_len_out.push_back(forward_len);
        reverse_pos_out.push_back(reverse_pos);
        reverse_len_out.push_back(reverse_len);
        return;
    }

    void add_unoriented(int pos1, const AlignData& data1, int pos2, const AlignData& data2, bool isfirst1) {
        if (!diagnostics) { return; }
        if (isfirst1) {
            ufirst_pos.push_back(pos1);
            ufirst_len.push_back(data1.len);
            usecond_pos.push_back(pos2);
            usecond_len.push_back(data2.len);
        } else {
            ufirst_pos.push_back(pos2);
            ufirst_len.push_back(data2.len);
            usecond_pos.push_back(pos1);
            usecond_len.push_back(data1.len);
        }
    }

    void add_onemapped(int pos, const AlignData& data) {
        if (!diagnostics) { return; }
        onemap_pos.push_back(pos);
        onemap_len.push_back(data.len);
        return;
    }

    void add_single(int pos, const AlignData& data) {
        if (!diagnostics) { return; }
        single_pos.push_back(pos);
        single_len.push_back(data.len);
        return;
    }

    void add_interchr(int pos, const AlignData& data, const char* name, bool isfirst) {
        if (!diagnostics) { return; }
        if (isfirst) { 
            ifirst_pos.push_back(pos);
            ifirst_len.push_back(data.len);
            interchr_names_1.push_back(std::string(name));
        } else {
            isecond_pos.push_back(pos);
            isecond_len.push_back(data.len);
            interchr_names_2.push_back(std::string(name));
        }
        return;
    }

    const bool diagnostics;
    int totals;
    int forward_pos, reverse_pos, forward_len, reverse_len;

    std::deque<int> forward_pos_out, forward_len_out, reverse_pos_out, reverse_len_out;
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

SEXP extract_pair_data(SEXP bam, SEXP index, SEXP chr, SEXP start, SEXP end, SEXP mapq, SEXP dedup, SEXP diagnostics) try {
    // Checking input values.
    if (!isInteger(mapq) || LENGTH(mapq)!=1) {
        throw std::runtime_error("mapping quality should be an integer scalar");
    }    
    const int minqual=asInteger(mapq);

    if (!isLogical(dedup) || LENGTH(dedup)!=1) {
        throw std::runtime_error("duplicate removal should be a logical scalar"); 
    }
    const bool rmdup=asLogical(dedup);

    if (!isLogical(diagnostics) || LENGTH(diagnostics)!=1) { 
        throw std::runtime_error("diagnostics specification should be a logical scalar"); 
    }
    const bool getnames=asLogical(diagnostics);

    // Initializing odds and ends.
    BamFile bf(bam, index);
    BamRead br;
    BamIterator biter(bf, chr, start, end);
    OutputContainer oc(getnames);
        
    typedef std::map<std::pair<int, std::string>, AlignData> Holder;
    std::deque<Holder> all_holders(4); // four holders, one for each strand/first combination; cut down searches.
    std::pair<int, std::string> current;
    Holder::iterator ith;
    int curpos, mate_pos;
    AlignData algn_data;
    bool am_mapped, is_first;

    bool mate_is_in;
    std::set<std::string> identical_pos;
    std::set<std::string>::iterator itip;
    int last_identipos=-1;

    while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){
        ++oc.totals;
        curpos = (br.read->core).pos + 1; // Getting 1-indexed position.
        br.extract_data(algn_data);
        am_mapped=br.is_well_mapped(minqual, rmdup);

        /* Reasons to not add a read: */
       
//        // If we can see that it is obviously unmapped (IMPOSSIBLE for a sorted file).
//        if (((br.read -> core).flag & BAM_FUNMAP)!=0) { 
//            // We don't filter by additional mapping criteria, as we need to search 'holder' to pop out the partner and to store diagnostics.
//            continue;
//        } 
        
        // If it's a singleton.
        if (((br.read -> core).flag & BAM_FPAIRED)==0) {
            if (am_mapped) { oc.add_single(curpos, algn_data); }
            continue;
        }

        // Or, if we can see that its partner is obviously unmapped.
        if (((br.read -> core).flag & BAM_FMUNMAP)!=0) {
            if (am_mapped) { oc.add_onemapped(curpos, algn_data); }
            continue;
        }

        // Or if it's inter-chromosomal.
        is_first=(((br.read->core).flag & BAM_FREAD1)!=0);
        if (is_first==(((br.read->core).flag & BAM_FREAD2)!=0)) { 
            std::stringstream err;
            err << "read '" << bam_get_qname(br.read) << "' must be either first or second in the pair";
            throw std::runtime_error(err.str()); 
        }
      
        if ((br.read -> core).mtid!=(br.read -> core).tid) { 
            if (am_mapped) { oc.add_interchr(curpos, algn_data, bam_get_qname(br.read), is_first); } 
            continue;
        }

        /* Checking the map and adding it if it doesn't exist. */
        
        current.second.assign(bam_get_qname(br.read));
        mate_pos = (br.read -> core).mpos + 1; // 1-indexed position, again.
        mate_is_in=false;
        if (mate_pos < curpos) {
            mate_is_in=true;
        } else if (mate_pos == curpos) {
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
            current.first = mate_pos;
            Holder& holder=all_holders[int(!is_first) + 2*int(bam_is_mrev(br.read))];
            ith=holder.find(current);

            if (ith != holder.end()) { 
                if (!am_mapped) {
                    // Searching to pop out the mate, to reduce the size of 'holder' for the remaining searches (and to store diagnostics).
                    oc.add_onemapped((ith->first).first, ith->second);
                    holder.erase(ith);
                    continue;
                }

                oc.add_genuine(curpos, algn_data, (ith->first).first, ith->second, is_first);
                holder.erase(ith);
            } else if (am_mapped) {
                // Only possible if the mate didn't get added because 'am_mapped' was false.
                oc.add_onemapped(curpos, algn_data);
            }
        } else if (am_mapped) {
            current.first = curpos;
            Holder& holder=all_holders[int(is_first) + 2*int(algn_data.is_reverse)];
            holder[current] = algn_data;
        }
    }

    // Leftovers treated as one_unmapped; marked as paired, but the mate is not in file.
    for (size_t h=0; h<all_holders.size(); ++h) { 
        Holder& holder=all_holders[h];
        for (ith=holder.begin(); ith!=holder.end(); ++ith) { 
            oc.add_onemapped((ith->first).first, ith->second);
        }
        holder.clear();
    }    

    // Storing all output.
    SEXP output=PROTECT(allocVector(VECSXP, getnames ? 9 : 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, 2));
        SEXP left=VECTOR_ELT(output, 0);
        store_int_output(left, 0, oc.forward_pos_out);
        store_int_output(left, 1, oc.forward_len_out);
        
        SET_VECTOR_ELT(output, 1, allocVector(VECSXP, 2));
        SEXP right=VECTOR_ELT(output, 1);
        store_int_output(right, 0, oc.reverse_pos_out);
        store_int_output(right, 1, oc.reverse_len_out);
    
        if (getnames) {
            SET_VECTOR_ELT(output, 2, ScalarInteger(oc.totals));
            
            SET_VECTOR_ELT(output, 3, allocVector(VECSXP, 2));
            SEXP singles=VECTOR_ELT(output, 3);
            store_int_output(singles, 0, oc.single_pos);
            store_int_output(singles, 1, oc.single_len);

            SET_VECTOR_ELT(output, 4, allocVector(VECSXP, 2));
            SEXP first=VECTOR_ELT(output, 4);
            store_int_output(first, 0, oc.ufirst_pos);
            store_int_output(first, 1, oc.ufirst_len);
            
            SET_VECTOR_ELT(output, 5, allocVector(VECSXP, 2));
            SEXP second=VECTOR_ELT(output, 5);
            store_int_output(second, 0, oc.usecond_pos);
            store_int_output(second, 1, oc.usecond_len);

            SET_VECTOR_ELT(output, 6, allocVector(VECSXP, 2));
            SEXP onemap=VECTOR_ELT(output, 6);
            store_int_output(onemap, 0, oc.onemap_pos);
            store_int_output(onemap, 1, oc.onemap_len);

            SET_VECTOR_ELT(output, 7, allocVector(VECSXP, 3));
            SEXP interchr1=VECTOR_ELT(output, 7);
            store_int_output(interchr1, 0, oc.ifirst_pos);
            store_int_output(interchr1, 1, oc.ifirst_len);
            store_names(interchr1, 2, oc.interchr_names_1);

            SET_VECTOR_ELT(output, 8, allocVector(VECSXP, 3));
            SEXP interchr2=VECTOR_ELT(output, 8);
            store_int_output(interchr2, 0, oc.isecond_pos);
            store_int_output(interchr2, 1, oc.isecond_len);
            store_names(interchr2, 2, oc.interchr_names_2);
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

/* Getting reads on other unprocessed chromosomes, unmapped reads. */

SEXP get_leftovers (SEXP bam, SEXP index, SEXP processed) try { 
    BamFile bf(bam, index);
    BamRead br;

    if (!isString(processed)) { throw std::runtime_error("names of processed chromosomes should be strings"); }
    const int nchr=LENGTH(processed);
    std::set<std::string> already_there;
    for (int i=0; i<nchr; ++i) {
        already_there.insert(std::string(CHAR(STRING_ELT(processed, i))));        
    }

    // Getting the reads mapped to chromosomes we didn't look at due to 'restrict'.
    int leftovers=0;
    std::set<std::string>::iterator iat;
    for (int cid=0; cid<bf.header->n_targets; ++cid) {
        iat=already_there.find(std::string(bf.header->target_name[cid]));
        if (iat!=already_there.end()) { continue; }
        BamIterator biter(bf, cid, 0, bf.header->target_len[cid]);
        while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){ ++leftovers; }
    } 
    
    // Also getting the unmapped guys. 
    BamIterator biter(bf);
    while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){ ++leftovers; }
    return(ScalarInteger(leftovers));
} catch (std::exception &e) {
    return mkString(e.what());
}

