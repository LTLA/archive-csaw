#include "csaw.h"

extern "C" {

/* We assume that incoming elements are sorted by chr -> start -> end. We then proceed 
 * to aggregate elements by chaining together elements that are less than 'tolerance'
 * apart and, if required, have the same 'sign'. We also split them if the difference
 * between the first start and the last end is greater than 'max_size'.
 */

SEXP R_merge(SEXP chrs, SEXP start, SEXP end, SEXP sign, SEXP tolerance, SEXP max_size) try {
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
	if (!n) { 
		SEXP output=PROTECT(NEW_LIST(4));
		try {
			for (int x=0; x<4; ++x) { SET_VECTOR_ELT(output, x, NEW_INTEGER(0)); }
		} catch (std::exception& e) {
			UNPROTECT(1);
			throw;
		}
		UNPROTECT(1);
		return output;
	}
		
	// Proceeding with the merge operation.
	SEXP output=PROTECT(NEW_LIST(4));
	try {
		SET_VECTOR_ELT(output, 0, NEW_INTEGER(n));
		int* optr=INTEGER_POINTER(VECTOR_ELT(output, 0));
		*optr=1;
		int last_end=*eptr, last_sign=*lptr, current_start=*sptr;
		for (int i=1; i<n; ++i) {
			optr[i]=optr[i-1];
			const bool diffchr=(cptr[i]!=cptr[i-1]);
			if (diffchr 											// Obviously, changing if we're on a different chromosome.
				|| sptr[i]-eptr[i-1]-1 > tol						// Space between windows, start anew if this is greater than the tolerance.
				|| lptr[i]!=lptr[i-1] 								// Checking if the sign is consistent.
				|| (limit_size && eptr[i]-current_start >= maxs) 	// Width is end-start+1, but '+1' gets absorbed when '>' turns into '>='.
		   	) { 
				++optr[i]; 
				current_start=sptr[i];
			}

			/* Note that any nested regions with the opposite sign as the parent will break
 			 * any stretch involving the parent. This is probably the correct interpretation
 			 * and it is also convenient to code. Mind you, these shouldn't really be observed
 			 * anyway, except maybe at the ends of the chromosome where trimming results in
 			 * everything having the same end point (and even then, you'd need a striped pattern
 			 * of changes to get alternating signs throughout). 
 			 */ 
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
			osptr[curgroup-1]=sptr[i];
			int& curend=(oeptr[curgroup-1]=eptr[i]);
			while (curgroup==optr[++i]) {
				if (curend < eptr[i]) { curend=eptr[i]; }
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
