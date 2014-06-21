#include "csaw.h"

#include <string>
#include <cstdio>
#include <map>

extern "C" {

/* This function spits out the ID for each exon, with some degree of strand-awareness, 
 * such that the first exon in the gene is labelled as exon 1, then 2, 3, etc. They
 * are assumed to be stored in genomic order so we just label them as-is.
 */

SEXP collate_exon_data (SEXP geneid, SEXP strand, SEXP start, SEXP end) try {
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
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what()); 
}

/* This function takes a 'start' and 'end' index. These two indices indicate the
 * range of overlaps for a particular annotation/region combination. The
 * annotation indices corresponding to each overlap are in 'indices', so
 * 'indices[start <= x < end]' will give integer indices that point to the gene
 * information i.e. 'symbols', 'features', 'strand'. The overlaps are also assigned
 * distances in 'dists', so 'dists[start <= x < end]' will give the distances from
 * the edge of the region to the annotated bit (if the pointer is not NULL).
 * The function will return a string collating all information for that annotated
 * feature (i.e., collate exon-level overlaps to a gene-level string).
 */

#define BUFSIZE 100
char make_int_box [BUFSIZE];
int errcode=0;
char * ssfake (const int& x) { // Don't want to load sstream because of Rf_length conflict.
//	if (x==NA_INTEGER || x==0) { // Had all these checks, but what kind of integer would have over 100 digits anyway? ~10 (+1 for the null, and another for the strand) would be the max.
//		;
//	} else if (x>0) { 
//		if (int(std::ceil(std::log10(double(x)))) + 2 > BUFSIZE) { // Checking for safety. 
//			throw std::runtime_error("insufficient space in the buffer for integer conversion"); 
//		}
//	} else {
//		throw std::runtime_error("negative integers should not present during annotation"); 
//	}
	errcode=std::sprintf(make_int_box, "%d", x); // Don't want to use snprintf, avoid C++11. 
	if (errcode >= 0) { return make_int_box; }
	throw std::runtime_error("error in string to integer conversion for annotation");
}

std::string digest2string (const int start, const int end, const int* indices, const int* dists, SEXP symbols, const int* features, const int* strand) {
	std::string ss(CHAR(STRING_ELT(symbols, indices[start])));
	ss += '|'; 

	// Deciding what to print.
	if (end==start) {
		throw std::runtime_error("empty vector of collected exon numbers"); 
	} else if (end==start+1) {
		if (features[indices[start]]==-1) {
			ss += "I"; 
		} else {
			ss += ssfake(features[indices[start]]);
		}
	} else {		
		int index=start;
		if (features[indices[index]]==-1) { ++index; } 
		ss += ssfake(features[indices[index]]);
		bool wasokay=false;

		// Running through and printing all stretches of contiguous exons.
		while ((++index)<end) {
			if (features[indices[index]]==features[indices[index-1]]+1) { 
				wasokay=true;
			} else {
				if (wasokay) {
					ss += '-';
 				    ss += ssfake(features[indices[index-1]]);
					wasokay=false;
				}
				ss += ',';
 			    ss += ssfake(features[indices[index]]);
			}
		}
		if (wasokay) {
 		    ss += '-';
			ss += ssfake(features[indices[end-1]]); 
		}
	}
	
	// Adding the strand and distance information.	
	ss += '|';
    ss += (strand[indices[start]] ? '+' : '-');
	if (dists!=NULL) {
		int lowest=dists[start];
		for (int index=start+1; index < end; ++index) {
			if (lowest > dists[index]) { lowest=dists[index]; }
		}
		ss += '[';
		ss += ssfake(lowest);
		ss += ']';
	}
	return ss;
}

SEXP annotate_overlaps (SEXP N, SEXP fullQ, SEXP fullS, SEXP leftQ, SEXP leftS, SEXP leftDist,
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
	int curend;
	int* curx_p;
	const int* curn_p;
	const int* cur_qptr;
	const int* cur_sptr;
	const int* cur_dptr;
	SEXP curout;

	for (int curreg=0; curreg<nin; ++curreg) {
		// Adding all overlaps of each type. Assuming that findOverlaps gives ordered output, which it should.
		std::string resultstr;
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
	
			// For the current region, we get everything in the current overlap; for each overlapping gene, we collate its features.
			int& curx=*curx_p;
			const int& curn=*curn_p;
			while (curx < curn && cur_qptr[curx]==curreg) {
				const int& curindex=cur_sptr[curx];
				if (curindex >= nsym) { throw std::runtime_error("symbol out of range for overlap index"); }
				curend=curx+1;
				while (curend < curn && cur_qptr[curend]==curreg && giptr[cur_sptr[curend]]==giptr[curindex]) { ++curend; }
				if (!resultstr.empty()) { resultstr += ","; }
				resultstr += digest2string(curx, curend, cur_sptr, cur_dptr, symbol, gfptr, gsptr);
				curx=curend;
			}
			SET_STRING_ELT(curout, curreg, mkChar(resultstr.c_str()));
			resultstr.clear();
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

}
