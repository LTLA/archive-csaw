#include "csaw.h"

extern "C" {

/* We assume that incoming elements are sorted by chr -> start -> end. We then proceed 
 * to aggregate elements by chaining together elements that are less than 'tolerance'
 * apart and, if required, have the same 'sign'. We also split them if the difference
 * between the first start and the last end is greater than 'max_size'.
 */

SEXP merge_windows(SEXP chrs, SEXP start, SEXP end, SEXP sign, SEXP tolerance, SEXP max_size) try {
	if (!IS_INTEGER(chrs)) { throw std::runtime_error("chromosomes should be a integer vector"); }
	if (!IS_INTEGER(start)) { throw std::runtime_error("start vector should be integer"); }
	if (!IS_INTEGER(end)) { throw std::runtime_error("end vector should be integer"); }
	if (!IS_LOGICAL(sign)) { throw std::runtime_error("sign vector should be logical"); }
	if (!IS_INTEGER(tolerance) || LENGTH(tolerance)!=1) { throw std::runtime_error("tolerance should be an integer scalar"); }

	// Setting up pointers.
	const int* cptr=INTEGER_POINTER(chrs);
	const int* sptr=INTEGER_POINTER(start);
	const int* eptr=INTEGER_POINTER(end);
	const int* lptr=LOGICAL_POINTER(sign);
	const int tol=INTEGER_VALUE(tolerance);
	
	// Checking whether we need to supply a maximum size.
	if (!IS_INTEGER(max_size) || LENGTH(max_size) > 1) { throw std::runtime_error("maximum size should be an integer scalar"); }
	const bool limit_size=(LENGTH(max_size)==1);
	const int maxs=(limit_size ? INTEGER_VALUE(max_size) : 0);
	
	// Providing some protection against an input empty list.
	const int n = LENGTH(chrs);
	if (n!=LENGTH(start) || n!=LENGTH(end) || n!=LENGTH(sign)) { throw std::runtime_error("lengths of vectors are not equal"); 	}
	if (n==0) { throw std::runtime_error("no elements provided for clustering"); }
		
	// Proceeding with the merge operation.
	SEXP output=PROTECT(NEW_LIST(4));
	try {
		SET_VECTOR_ELT(output, 0, NEW_INTEGER(n));
		int* optr=INTEGER_POINTER(VECTOR_ELT(output, 0));
		*optr=1;
		int current_start=*sptr, last_end=*eptr;
		bool diffchr, diffsign;

		for (int i=1; i<n; ++i) {
			optr[i]=optr[i-1];
			diffchr=(cptr[i]!=cptr[i-1]);
 		   	diffsign=(lptr[i]!=lptr[i-1]);

			if (diffchr 											// Obviously, changing if we're on a different chromosome.
				|| sptr[i]-last_end-1 > tol							// Space between windows, start anew if this is greater than the tolerance.
				|| diffsign 										// Checking if the sign is consistent.
				|| (limit_size && eptr[i]-current_start >= maxs) 	// Width is end-start+1, but '+1' gets absorbed when '>' turns into '>='.
		   	) { 
				++optr[i]; 
				current_start=sptr[i];
			}

			/* Fully nested regions don't have a properly defined interpretation when it comes
 			 * to splitting things by sign. We only support consideration of nested regions where
 			 * either of the boundaries are the same. That can be considered to break the 
 			 * previous stretch if it had opposite sign. Otherwise, the next window would have to
 			 * make the decision to match the parent window or its nested child.
 			 *
 			 * If the start is the same, the window with the earlier end point should be ordered
 			 * first, so that won't pass the first 'if'. This means that it'll only enter with
 			 * a fully nested window. Start and end-point equality might be possible at the ends 
 			 * of chromosomes where trimming enforces sameness, but full nesting should not be observed.
 			 *
 			 * If the nested region has the same sign as the parent, then everything proceeds
 			 * normally i.e. same cluster. We make sure to keep 'last_end' as the parent end in
 			 * such cases. This ensures that the following windows get a change to compute
 			 * distances to the parent end (which should be closer).
 			 */
 		    if (!diffchr && eptr[i] < eptr[i-1]) {
 			   	if (diffsign) { throw std::runtime_error("fully nested windows of opposite sign are not supported"); } 
			} else { last_end=eptr[i]; }
		}

		// Now, identifying the chromosome, start and end of each region.
		const int ngroups=optr[n-1];
		SET_VECTOR_ELT(output, 1, NEW_INTEGER(ngroups));
		int* ocptr=INTEGER_POINTER(VECTOR_ELT(output, 1));
		SET_VECTOR_ELT(output, 2, NEW_INTEGER(ngroups));
		int* osptr=INTEGER_POINTER(VECTOR_ELT(output, 2));
		SET_VECTOR_ELT(output, 3, NEW_INTEGER(ngroups));
		int* oeptr=INTEGER_POINTER(VECTOR_ELT(output, 3));
		int i=0; 
		while (i<n) {
			const int& curgroup=optr[i];
			ocptr[curgroup-1]=cptr[i];
			osptr[curgroup-1]=sptr[i]; // Sorted by sign, remember.
			int& curend=(oeptr[curgroup-1]=eptr[i]);
			++i;
			while (i < n && curgroup==optr[i]) {
				if (curend < eptr[i]) { curend=eptr[i]; }
				++i;
			}
		}
	} catch (std::exception& e){
		UNPROTECT(1);
		throw;
	}
	
	UNPROTECT(1);
	return output;		
} catch (std::exception& e) {
	return mkString(e.what());
}

}
