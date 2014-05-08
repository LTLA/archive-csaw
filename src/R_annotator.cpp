#include "csaw.h"
#include <string>
#include <sstream>
#include <map>

extern "C" {

/* This function spits out the ID for each exon, with some degree of strand-awareness, 
 * such that the first exon in the gene is labelled as exon 1, then 2, 3, etc. They
 * are assumed to be stored in genomic order so we just label them as-is.
 */

SEXP R_collate_exon_data (SEXP geneid, SEXP strand, SEXP start, SEXP end) try {
	// Checking inputs.
	if (!IS_INTEGER(geneid)) { throw std::runtime_error("gene ID vector should be integer"); }
	if (!IS_LOGICAL(strand)) { throw std::runtime_error("vector of strands should be logical"); }
	if (!IS_INTEGER(start) || !IS_INTEGER(end)) { throw std::runtime_error("start/end positions and indices should be integer vectors"); }
	const int n=LENGTH(geneid);
	if (n!=LENGTH(strand)) { throw std::runtime_error("strand/ID vectors should have same length"); }
	if (n!=LENGTH(start) || n!=LENGTH(end)) { throw std::runtime_error("start/end/index vectors should have the same length"); }
	const int* gixptr=INTEGER_POINTER(geneid),
		* strptr=LOGICAL_POINTER(strand),
		* staptr=INTEGER_POINTER(start),
		* endptr=INTEGER_POINTER(end);
	sort_row_index<int> endcomp(endptr);
	
	// Scanning through to determine the number of unique genes.
	int nuniq=1;
	for (int x=1; x<n; ++x) {
		if (gixptr[x]!=gixptr[x-1]) { 
			++nuniq; 
		} else if (strptr[x]!=strptr[x-1]) { 
			throw std::runtime_error("exons of the same gene should have the same strand"); 
		} else if (staptr[x]<staptr[x-1]) {
			throw std::runtime_error("exons of the same gene should be sorted by the start index"); 
		}	
	}

	// Setting up output structures.
	SEXP output=PROTECT(NEW_LIST(2));
try {
	SET_VECTOR_ELT(output, 0, NEW_INTEGER(n));
	int* eiptr=INTEGER_POINTER(VECTOR_ELT(output, 0));
	SET_VECTOR_ELT(output, 1, NEW_LIST(3));
	SEXP genebody=VECTOR_ELT(output, 1);
	SET_VECTOR_ELT(genebody, 0, NEW_INTEGER(nuniq));
	SET_VECTOR_ELT(genebody, 1, NEW_INTEGER(nuniq));
	SET_VECTOR_ELT(genebody, 2, NEW_INTEGER(nuniq));
	int * oiptr=INTEGER_POINTER(VECTOR_ELT(genebody, 0)),
		* osptr=INTEGER_POINTER(VECTOR_ELT(genebody, 1)),
		* oeptr=INTEGER_POINTER(VECTOR_ELT(genebody, 2));

	int curex=0, counter=0, ngene=0;
	std::deque<int> current_indices;
	while (curex<n) {
		const int& current=gixptr[curex];
		oiptr[ngene]=curex+1;
		osptr[ngene]=staptr[curex]; // Input should be sorted by start positions within each gene ID, so first element should be earliest.
		int& last_end=(oeptr[ngene]=endptr[curex]);
		++ngene;

		// Resorting by end positions, if the gene is on the negative strand.
		if (!strptr[curex]) { 
			do {
				current_indices.push_back(curex);
				++curex;
			} while (curex < n && current==gixptr[curex]);
			
			const int& stretch=current_indices.size();
			std::sort(current_indices.begin(), current_indices.end(), endcomp); 
			for (counter=0; counter<stretch; ++counter) { 
				eiptr[current_indices[counter]]=stretch - counter;
			}
			last_end=endptr[current_indices.back()];
			current_indices.clear();
		} else {
			counter=1;
			do {
				eiptr[curex]=counter;
				if (last_end < endptr[curex]) { last_end=endptr[curex]; }
				++counter;
				++curex;
			} while (curex < n && current==gixptr[curex]);
		}
	}
} catch (std::exception& e) {
	UNPROTECT(1);
	throw e;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what()); 
}

/* This function takes a region and finds all features that start or end in
 * that region. To wit; the start and end is re-defined according to the
 * strand, so it will process them correctly. If both the start and the end
 * match, it'll only keep the closer value.
 */

struct overlaps { 
	overlaps (int o, int a, int b, int c=-1) : original(o), min(a), max(b), dist(c) {}
	int original, min, max, dist; 
};
std::string digest2string (std::map<int, overlaps>& tango, SEXP symbols, const int* gsptr) {
	std::string temp;
	for (std::map<int, overlaps>::iterator it=tango.begin(); it!=tango.end(); ++it) {
		std::stringstream ss;
		const int& curid = (it->second).original;
		ss << CHAR(STRING_ELT(symbols, curid)) << '|'; 

		// Deciding what to print.
		if ((it->second).min >=0) { 
			ss << (it->second).min;
			if ((it->second).min!=(it->second).max) { ss << '-' << (it->second).max; }
		} else {
			ss << "I";
		}
		ss << '|' << (gsptr[curid] ? '+' : '-');
		if ((it->second).dist > 0) { ss << "[" << (it->second).dist << "]"; }
		if (!temp.empty()) { temp += ","; }
		temp += ss.str();
	}
	return temp;
}

SEXP R_annotate (SEXP N, SEXP fullQ, SEXP fullS, SEXP leftQ, SEXP leftS, SEXP leftDist,
		SEXP rightQ, SEXP rightS, SEXP rightDist, 
		SEXP symbol, SEXP genefeature, SEXP geneid, SEXP genestr) try {
	if (!IS_INTEGER(N) || LENGTH(N)!=1) { throw std::runtime_error("N should be a integer scalar"); }
	const int nin=INTEGER_VALUE(N);
	if (!IS_INTEGER(fullQ) || !IS_INTEGER(fullS)) { throw std::runtime_error("full overlap query/subject IDs should be integer vectors"); }
	const int nfull=LENGTH(fullQ);
	if (nfull!=LENGTH(fullS)) { throw std::runtime_error("full overlap vectors should have equal length"); }
	if (!IS_INTEGER(leftQ) || !IS_INTEGER(leftS) || !IS_INTEGER(leftDist)) { throw std::runtime_error("left overlap query/subject/distances should be integer vectors"); }
	const int nleft=LENGTH(leftQ);
	if (nleft!=LENGTH(leftS) || nleft!=LENGTH(leftDist)) { throw std::runtime_error("left overlap vectors should have equal length"); }
	if (!IS_INTEGER(rightQ) || !IS_INTEGER(rightS) || !IS_INTEGER(rightDist)) { throw std::runtime_error("right overlap query/subject/distances should be integer vectors"); }
	const int nright=LENGTH(rightQ);
	if (nright!=LENGTH(rightS) || nright!=LENGTH(rightDist)) { throw std::runtime_error("right overlap vectors should have equal length"); }

	// Declaring metafeatures.
	if (!IS_CHARACTER(symbol)) { throw std::runtime_error("symbols should be a character vector"); }
	const int nsym=LENGTH(symbol);
	if (!IS_INTEGER(geneid)) { throw std::runtime_error("gene IDs should be an integer vector"); }
	if (!IS_INTEGER(genefeature)) { throw std::runtime_error("gene feature should be an integer vector"); }
	if (!IS_LOGICAL(genestr)) {throw std::runtime_error("gene strand should be a logical vector"); }
	if (nsym!=LENGTH(geneid) || nsym!=LENGTH(genefeature) || nsym!=LENGTH(genestr)) { 
		throw std::runtime_error("gene data vectors should have the same length"); 
	}

	// Setting up pointers.
	const int * fqptr=INTEGER_POINTER(fullQ),
		  * fsptr=INTEGER_POINTER(fullS),
		  * lqptr=INTEGER_POINTER(leftQ),
		  * lsptr=INTEGER_POINTER(leftS),
		  * ldptr=INTEGER_POINTER(leftDist),
		  * rqptr=INTEGER_POINTER(rightQ),
		  * rsptr=INTEGER_POINTER(rightS),
		  * rdptr=INTEGER_POINTER(rightDist),
		  * giptr=INTEGER_POINTER(geneid),
		  * gfptr=INTEGER_POINTER(genefeature),
		  * gsptr=LOGICAL_POINTER(genestr);

	// Okay, now going through them and assembling the output vectors. 
	SEXP output=PROTECT(NEW_LIST(3));
try {
	SET_VECTOR_ELT(output, 0, NEW_CHARACTER(nin));
	SEXP full_out=VECTOR_ELT(output, 0);
	SET_VECTOR_ELT(output, 1, NEW_CHARACTER(nin));
	SEXP left_out=VECTOR_ELT(output, 1);
	SET_VECTOR_ELT(output, 2, NEW_CHARACTER(nin));
	SEXP right_out=VECTOR_ELT(output, 2);

	int fullx=0, leftx=0, rightx=0;
	std::map<int, overlaps> processed;
	std::map<int, overlaps>::iterator itp;
	for (int curreg=0; curreg<nin; ++curreg) {

		// Adding the full overlaps. Assuming that findOverlaps gives ordered output, which it should.
		while (fullx < nfull && fqptr[fullx]==curreg) {
			const int& curindex=fsptr[fullx];
			if (curindex >= nsym) { throw std::runtime_error("symbol index out of range for full overlaps"); }
					
			// Smart insert. A gene body overlap is the baseline (-1 for curfeature), but anything else replaces it.
			const int& curgene=giptr[curindex];
			const int& curfeature=gfptr[curindex];
			itp=processed.lower_bound(curgene);
			if (itp==processed.end() || processed.key_comp()(curgene, itp->first)) {
				processed.insert(itp, std::make_pair(curgene, overlaps(curindex, curfeature, curfeature)));
			} else {
				if (curfeature >= 0) { 
					if (curfeature < (itp->second).min || (itp->second).min<0) { (itp->second).min = curfeature; }
					if (curfeature > (itp->second).max || (itp->second).max<0) { (itp->second).max = curfeature; }
				}
			}
			++fullx;
		}

		std::string temp=digest2string(processed, symbol, gsptr);
//		std::cout << curreg << "\t" << temp << std::endl;
		SET_STRING_ELT(full_out, curreg, mkChar((digest2string(processed, symbol, gsptr)).c_str()));
//		std::cout << "Done!" << std::endl;
		processed.clear();

		// Adding the left and right overlaps. The key thing here is that only the closest feature is kept.
		for (int mode=0; mode<2; ++mode) { 
			int& curx = (mode ? leftx : rightx);
			const int & ncur = (mode ? nleft : nright);
			const int* cur_qptr = (mode ? lqptr : rqptr);
			const int* cur_sptr = (mode ? lsptr : rsptr);
			const int* cur_dptr = (mode ? ldptr : rdptr);

			while (curx < ncur && cur_qptr[curx]==curreg) {
				const int& curdist=cur_dptr[curx];
				if (curdist > 0) { 
					const int& curindex=cur_sptr[curx];
					if (curindex >= nsym) { throw std::runtime_error("symbol index out of range for left/right overlaps"); }
		
					// Smart insert. We don't bother with intronic overlaps on the left and right.
					const int& curgene=giptr[curindex];
					const int& curfeature=gfptr[curindex];
					if (curfeature >= 0)  {
						itp=processed.lower_bound(curgene);
						if (itp==processed.end() || processed.key_comp()(curgene, itp->first)) {
							processed.insert(itp, std::make_pair(curgene, overlaps(curindex, curfeature, curfeature, curdist)));
						} else {
							if (curfeature < (itp->second).min) { (itp->second).min = curfeature; }
							if (curfeature > (itp->second).max) { (itp->second).max = curfeature; }
							if (curfeature < (itp->second).dist) { (itp->second).dist = curdist; }
						}
					}
				}
				++curx;
			}

//			std::cout << (mode ? "LEFT:" : "RIGHT:") << "\t" << temp << std::endl;
			SET_STRING_ELT((mode ? left_out : right_out), curreg, mkChar((digest2string(processed, symbol, gsptr)).c_str()));
//			std::cout << "Done!" << std::endl;
			processed.clear();
		}
	}
} catch (std::exception& e) {
	UNPROTECT(1);
	throw e;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

}
