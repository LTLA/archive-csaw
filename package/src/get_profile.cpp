#include "csaw.h"

#ifdef DEBUG
#include <map>
#endif

/* This function scans through the track and pulls out local maxima. */

SEXP get_profile(SEXP starts, SEXP ends, SEXP regstarts, SEXP weights, SEXP range) try {

	if (!isInteger(starts) || !isInteger(ends)) { throw std::runtime_error("fragment start/end positions should be integer"); }
	const int nfrags=LENGTH(starts);
	if (LENGTH(ends)!=nfrags) { throw std::runtime_error("fragment start/end vectors should have same length"); }
	if (!isInteger(regstarts)) { throw std::runtime_error("region start/end positions should be integer"); }
	const int nregs=LENGTH(regstarts);
	if (nregs==0) { throw std::runtime_error("no regions supplied"); }
	if (!isReal(weights))  { throw std::runtime_error("weight vector should be double-precision"); }
	if (LENGTH(weights)!=nregs) { throw std::runtime_error("weight vector should have length equal to number of regions"); }

	// Setting up constructs.
	if (!isInteger(range) || LENGTH(range)!=1) { throw std::runtime_error("range distance should be an integer scalar"); }
	const int maxrange=asInteger(range);
	const int* fsptr=INTEGER(starts),
		  *feptr=INTEGER(ends),
		  *rsptr=INTEGER(regstarts);
	const double* wptr=REAL(weights);

	/* Setting up a separate profile for each region. This is necessary
	 * to ensure that the calculations are integer (despite weighting),
	 * in order to preserve numerical stability.
	 */
	const int totallen=2*maxrange+1;
	std::deque<int*> all_profiles(nregs);
	for (int curreg=0; curreg<nregs; ++curreg) {
		all_profiles[curreg]=(int*)R_alloc(totallen, sizeof(int));
		all_profiles[curreg]+=maxrange; // so index of 0 = distance of 0.
	}

	/* Running through the reads. We use a strategy of identifying the
	 * regions for each read and adding that to the profile, rather than
	 * setting up some queue (which would involve lots of insertions/deletions).
	 */
	const int* ptr=NULL, *ptr_copy=NULL, *terminal_s=rsptr+nregs;
	int dist, dist2;
	int* curprof;
	for (int frag=0; frag<nfrags; ++frag) {
		const int& curstart=fsptr[frag];
		const int& curend=feptr[frag];

		/* Getting all regions starting after the fragment. Don't bother
		 * trying to optimize with reducing the binary search by sorting
		 * the fragments beforehand; the earlier sort would take more time
		 * anyway, because it's nfrag*log(nfrag), not nfrag*log(nregs).
		 */
		ptr_copy=ptr=std::upper_bound(rsptr, terminal_s, curend);
		while (ptr!=terminal_s) {
			dist = *ptr - curend;
			if (dist > maxrange) { break; }
			curprof=all_profiles[ptr-rsptr];
			--(curprof[-dist+1]);
			dist2 = *ptr - curstart;
			++(curprof[dist2 >= maxrange ? -maxrange : -dist2]);
			++ptr;
		}

		// Getting all regions starting before the fragment end.
		ptr=ptr_copy;
		while (ptr!=rsptr) {
			--ptr;
			dist = curstart - *ptr;
			if (dist > maxrange) { break; }
			curprof=all_profiles[ptr-rsptr];
			++(curprof[dist < -maxrange ? -maxrange : dist]);
			dist2 = curend - *ptr;
			if (dist2 < maxrange) { --(curprof[dist2+1]); }
		}
	}

	// Setting up output.
	SEXP output=PROTECT(allocVector(REALSXP, totallen));
try{
	double* optr=REAL(output);
	for (int i=0; i<totallen; ++i) { optr[i]=0; }
	optr += maxrange; // 0 is now distance of zero.

	for (int curreg=0; curreg<nregs; ++curreg) {
		curprof=all_profiles[curreg];
		const double& curweight=wptr[curreg];
		optr[-maxrange]+=curprof[-maxrange]*curweight;
		for (int i=-maxrange+1; i<maxrange; ++i) {
			curprof[i]+=curprof[i-1]; // Compiling profile based on addition/subtraction instructions
			optr[i]+=curprof[i]*curweight;
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

