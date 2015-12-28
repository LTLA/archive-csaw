#include "csaw.h"

SEXP get_rle_counts(SEXP start, SEXP end, SEXP nr, SEXP space, SEXP first) try {
	if (!isInteger(nr) || LENGTH(nr)!=1) {  throw std::runtime_error("number of rows must be an integer scalar"); }
	if (!isInteger(space) || LENGTH(space)!=1) { throw std::runtime_error("spacing must be an integer scalar"); }
	if (!isLogical(first) || LENGTH(first)!=1) { throw std::runtime_error("decision to use first point must be a logical scalar"); }
	if (!isInteger(start))  { throw std::runtime_error("start vector must be integer"); }
	if (!isInteger(end))  { throw std::runtime_error("start vector must be integer"); }

	const int n=LENGTH(start);
	if (n!=LENGTH(end)) { throw std::runtime_error("start/end vectors must have equal length"); }
	const int nrows=asInteger(nr),
	  usefirst=asLogical(first),
	  spacing=asInteger(space);
	const int* sptr=INTEGER(start);
	const int* eptr=INTEGER(end);	

	SEXP output=PROTECT(allocVector(INTSXP, nrows));
	if (nrows==0) { 
		UNPROTECT(1);
		return output;
	}
	try {
		int* optr=INTEGER(output);
		for (int i=0; i<nrows; ++i) { optr[i]=0; }
		int left, right;
		for (int i=0; i<n; ++i) {
			// Get the zero-index corresponding to the smallest spacing point larger than the current inclusive start/end.
			if (eptr[i] < sptr[i])  { throw std::runtime_error("invalid coordinates for read start/ends"); }
			left=(sptr[i] < 2 ? 0 : int((sptr[i]-2)/spacing)+usefirst);
			right=(eptr[i] < 1 ? 0 : int((eptr[i]-1)/spacing)+usefirst);
			
			// Adding the steps for addition (at the start) and deletion (after the end) of each read.
			if (left<right) { 
				if (left<nrows) { ++optr[left]; }
				if (right<nrows) { --optr[right]; }
			}
		}

		// Running and computing the RLE, given the steps at each position.
		int cum=0;
		for (int i=0; i<nrows; ++i) { 
			cum+=optr[i];
			optr[i]=cum;
		}
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}
	UNPROTECT(1);
	return output;
} catch (std::exception& e){
	return mkString(e.what());
}
