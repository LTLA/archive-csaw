#include "csaw.h"

extern "C" {

SEXP R_get_rle_counts(SEXP start, SEXP end, SEXP nr, SEXP space) try {
	if (!IS_INTEGER(nr) || LENGTH(nr)!=1) {  throw std::runtime_error("number of rows must be an integer scalar"); }
	if (!IS_INTEGER(space) || LENGTH(space)!=1) { throw std::runtime_error("spacing must be an integer scalar"); }
	if (!IS_INTEGER(start))  { throw std::runtime_error("start vector must be integer"); }
	if (!IS_INTEGER(end))  { throw std::runtime_error("start vector must be integer"); }

	const int n=LENGTH(start);
	if (n!=LENGTH(end)) { throw std::runtime_error("start/end vectors must have equal length"); }
	const int nrows=INTEGER_VALUE(nr);
	const int spacing=INTEGER_VALUE(space);
	const int* sptr=INTEGER_POINTER(start);
	const int* eptr=INTEGER_POINTER(end);	

	SEXP output=PROTECT(NEW_INTEGER(nrows));
	try {
		int* optr=INTEGER_POINTER(output);
		for (int i=0; i<nrows; ++i) { optr[i]=0; }
		for (int i=0; i<n; ++i) {
			// Get the zero-index corresponding to the smallest spacing point larger than the current inclusive start/end.
			if (sptr[i] < 1 || eptr[i] < sptr[i])  { throw std::runtime_error("invalid coordinates for read start/ends"); }
			const int left=(sptr[i] < 2 ? 0 : int((sptr[i]-2)/spacing)+1);
			const int right=int((eptr[i]-1)/spacing)+1;

			// Adding the steps for addition (at the start) and deletion (after the end) of each read.
			if (left<right) { 
				++optr[left]; 
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

}
