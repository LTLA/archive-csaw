#include <sstream> // Before R.h in csaw.h, to resolve remapping issues (see Section 6, Writing R extensions).
#include "csaw.h"
#include <string>
#include <map>

/* This function spits out the ID for each exon, with some degree of strand-awareness, 
 * such that the first exon in the gene is labelled as exon 1, then 2, 3, etc. They
 * are assumed to be stored in genomic order so we just label them as-is.
 */

SEXP collate_exon_data (SEXP geneid, SEXP strand, SEXP start, SEXP end) try {
	// Checking inputs.
	if (!isInteger(geneid)) { throw std::runtime_error("gene ID vector should be integer"); }
	if (!isLogical(strand)) { throw std::runtime_error("vector of strands should be logical"); }
	if (!isInteger(start) || !isInteger(end)) { throw std::runtime_error("start/end positions and indices should be integer vectors"); }
	const int n=LENGTH(geneid);
	if (n!=LENGTH(strand)) { throw std::runtime_error("strand/ID vectors should have same length"); }
	if (n!=LENGTH(start) || n!=LENGTH(end)) { throw std::runtime_error("start/end/index vectors should have the same length"); }
	const int* gixptr=INTEGER(geneid),
		* strptr=LOGICAL(strand),
		* staptr=INTEGER(start),
		* endptr=INTEGER(end);
	sort_row_index<int> endcomp(endptr);
	
	// Scanning through to determine the number of unique genes.
	if (n==0) { throw std::runtime_error("no genes supplied for exonic aggregation"); }
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
	SEXP output=PROTECT(allocVector(VECSXP, 2));
try {
	SET_VECTOR_ELT(output, 0, allocVector(INTSXP, n));
	int* eiptr=INTEGER(VECTOR_ELT(output, 0));
	SET_VECTOR_ELT(output, 1, allocVector(VECSXP, 3));
	SEXP genebody=VECTOR_ELT(output, 1);
	SET_VECTOR_ELT(genebody, 0, allocVector(INTSXP, nuniq));
	SET_VECTOR_ELT(genebody, 1, allocVector(INTSXP, nuniq));
	SET_VECTOR_ELT(genebody, 2, allocVector(INTSXP, nuniq));
	int * oiptr=INTEGER(VECTOR_ELT(genebody, 0)),
		* osptr=INTEGER(VECTOR_ELT(genebody, 1)),
		* oeptr=INTEGER(VECTOR_ELT(genebody, 2));

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
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what()); 
}

/* This function collapses indices into a string. The overlaps are also assigned
 * distances in 'dists', so 'dists[start <= x < end]' will give the distances from
 * the edge of the region to the annotated bit (if the pointer is not NULL).
 * The function will return a string collating all information for that annotated
 * feature (i.e., collate exon-level overlaps to a gene-level string).
 */

std::string digest2string (const std::deque<int>& indices, const int* geneid, SEXP symbols, const int* features, const int* strand, const int* dists) {
	if (!indices.size()) { return ""; }
	std::stringstream ss;
	size_t start=0, end, index=0;

	while (start < indices.size()) {
		if (start!=0) { ss << ","; }
		ss << CHAR(STRING_ELT(symbols, indices[start])) << '|'; 
		end=start+1;	
		while (end < indices.size() && 	geneid[indices[end]]==geneid[indices[start]]) { ++end; }

		// Deciding what to print.
		if (end==start+1) {
			if (features[indices[start]]==-1) {
				ss << "I"; 
			} else {
				ss << features[indices[start]];
			}
		} else {	
			index=start;
			if (features[indices[start]]==-1) { ++index; }
			ss << features[indices[index]];
			bool wasokay=false;

			// Running through and printing all stretches of contiguous exons.
			while ((++index) < end) {
				if (features[indices[index]]==features[indices[index-1]]+1) { 
					wasokay=true;
				} else {
					if (wasokay) {
						ss << '-' << features[indices[index-1]];
						wasokay=false;
					}
					ss << ',' << features[indices[index]];
				}
			}
			if (wasokay) { ss << '-' << features[indices[index-1]]; }
		}
		
		// Adding the strand and distance information.	
		ss << '|' << (strand[indices[start]] ? '+' : '-');
		if (dists!=NULL) {
			int lowest=dists[indices[start]];
			for (index=start+1; index < end; ++index) {
				if (lowest > dists[indices[index]]) { lowest=dists[indices[index]]; }
			}
			ss << '[' << lowest << ']';
		}
		start=end;
	}
	return ss.str();
}

/* A sorting class for two elements. */

struct sort_pair_int_index { 
	sort_pair_int_index(const int* p1, const int* p2) : ptr1(p1), ptr2(p2) {}
	bool operator() (const int& l, const int& r) const { 
		if (ptr1[l]==ptr1[r]) { 
			if (ptr2[l]==ptr2[r]) { 
				return (l < r); 
			} else {
				return (ptr2[l] < ptr2[r]);
			}
		} else { 
			return (ptr1[l] < ptr1[r]); 
		}
	}
private:
	const int* ptr1, *ptr2;
};

/* The main function */

SEXP annotate_overlaps (SEXP N, SEXP fullQ, SEXP fullS, SEXP leftQ, SEXP leftS, SEXP leftDist,
		SEXP rightQ, SEXP rightS, SEXP rightDist, 
		SEXP symbol, SEXP genefeature, SEXP geneid, SEXP genestr) try {
	if (!isInteger(N) || LENGTH(N)!=1) { throw std::runtime_error("N should be a integer scalar"); }
	const int nin=asInteger(N);
	if (!isInteger(fullQ) || !isInteger(fullS)) { throw std::runtime_error("full overlap query/subject IDs should be integer vectors"); }
	const int nfull=LENGTH(fullQ);
	if (nfull!=LENGTH(fullS)) { throw std::runtime_error("full overlap vectors should have equal length"); }
	if (!isInteger(leftQ) || !isInteger(leftS) || !isInteger(leftDist)) { throw std::runtime_error("left overlap query/subject/distances should be integer vectors"); }
	const int nleft=LENGTH(leftQ);
	if (nleft!=LENGTH(leftS) || nleft!=LENGTH(leftDist)) { throw std::runtime_error("left overlap vectors should have equal length"); }
	if (!isInteger(rightQ) || !isInteger(rightS) || !isInteger(rightDist)) { throw std::runtime_error("right overlap query/subject/distances should be integer vectors"); }
	const int nright=LENGTH(rightQ);
	if (nright!=LENGTH(rightS) || nright!=LENGTH(rightDist)) { throw std::runtime_error("right overlap vectors should have equal length"); }

	// Declaring metafeatures.
	if (!isString(symbol)) { throw std::runtime_error("symbols should be a character vector"); }
	const int nsym=LENGTH(symbol);
	if (!isInteger(geneid)) { throw std::runtime_error("gene IDs should be an integer vector"); }
	if (!isInteger(genefeature)) { throw std::runtime_error("gene feature should be an integer vector"); }
	if (!isLogical(genestr)) {throw std::runtime_error("gene strand should be a logical vector"); }
	if (nsym!=LENGTH(geneid) || nsym!=LENGTH(genefeature) || nsym!=LENGTH(genestr)) { 
		throw std::runtime_error("gene data vectors should have the same length"); 
	}

	// Setting up pointers.
	const int * fqptr=INTEGER(fullQ),
		  * fsptr=INTEGER(fullS),
		  * lqptr=INTEGER(leftQ),
		  * lsptr=INTEGER(leftS),
		  * ldptr=INTEGER(leftDist),
		  * rqptr=INTEGER(rightQ),
		  * rsptr=INTEGER(rightS),
		  * rdptr=INTEGER(rightDist),
		  * giptr=INTEGER(geneid),
		  * gfptr=INTEGER(genefeature),
		  * gsptr=LOGICAL(genestr);
	
	// Okay, now going through them and assembling the output vectors. 
	SEXP output=PROTECT(allocVector(VECSXP, 3));
try {
	SET_VECTOR_ELT(output, 0, allocVector(STRSXP, nin));
	SEXP full_out=VECTOR_ELT(output, 0);
	SET_VECTOR_ELT(output, 1, allocVector(STRSXP, nin));
	SEXP left_out=VECTOR_ELT(output, 1);
	SET_VECTOR_ELT(output, 2, allocVector(STRSXP, nin));
	SEXP right_out=VECTOR_ELT(output, 2);

	int fullx=0, leftx=0, rightx=0;
	int* curx_p;
	const int* curn_p;
	const int* cur_qptr;
	const int* cur_sptr;
	const int* cur_dptr;
	SEXP curout;

	std::deque<int> allindices;
	sort_pair_int_index indexcomp(giptr, gfptr);
	int* distance_holder=(int*)R_alloc(nsym, sizeof(int));
	std::string resultstr;

	for (int curreg=0; curreg<nin; ++curreg) {
		// Adding all overlaps of each type. Assuming that findOverlaps gives ordered output, which it should.
		for (int mode=0; mode<3; ++mode) {
			if(mode==0) {
				curx_p=&fullx;
				curn_p=&nfull;
				cur_qptr=fqptr;
				cur_sptr=fsptr;
				cur_dptr=NULL;
				curout=full_out;
			} else if (mode==1) {
				curx_p=&leftx;
				curn_p=&nleft;
				cur_qptr=lqptr;
				cur_sptr=lsptr;
				cur_dptr=ldptr;
				curout=left_out;
			} else if (mode==2) {
				curx_p=&rightx;
				curn_p=&nright;
				cur_qptr=rqptr;
				cur_sptr=rsptr;
				cur_dptr=rdptr;
				curout=right_out;
			}

			// For the current region, we get everything in the current overlap.
			int& curx=*curx_p;
			const int& curn=*curn_p;
			allindices.clear();
			while (curx < curn && cur_qptr[curx]==curreg) {
				allindices.push_back(cur_sptr[curx]);
				if (allindices.back() >= nsym) { throw std::runtime_error("symbol out of range for overlap index"); }
				if (cur_dptr!=NULL) { distance_holder[allindices.back()] = cur_dptr[curx]; } // Storing distance to each overlapped feature.
				++curx;
			}

			// Sorting by gene index, then feature index; then collapsing into a string.
			std::sort(allindices.begin(), allindices.end(), indexcomp);
			resultstr = digest2string(allindices, giptr, symbol, gfptr, gsptr, (cur_dptr!=NULL ? distance_holder : NULL));
			SET_STRING_ELT(curout, curreg, mkChar(resultstr.c_str()));
		} 
	}
} catch (std::exception& e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

