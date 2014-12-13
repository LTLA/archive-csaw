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

	// Setting up output.
	const int totallen=2*maxrange+1;
	SEXP output=PROTECT(allocVector(REALSXP, totallen));
try{ 
	double* optr=REAL(output);
	for (int i=0; i<totallen; ++i) { optr[i]=0; }
	optr += maxrange; // 0 is now distance of zero.

	// Running through the reads.
	const int* ptr=NULL, *ptr_copy=NULL, *terminal_s=rsptr+nregs;
	int dist, dist2;
	for (int frag=0; frag<nfrags; ++frag) {
		const int& curstart=fsptr[frag];
		const int& curend=feptr[frag];

		// Getting all regions starting after the fragment.
		ptr_copy=ptr=std::upper_bound(rsptr, terminal_s, curend);
		while (ptr!=terminal_s) { 
			dist = *ptr - curend;
			if (dist > maxrange) { break; }
			optr[-dist+1] -= wptr[ptr-rsptr];
			dist2 = *ptr - curstart;
			optr[dist2 >= maxrange ? -maxrange : -dist2] += wptr[ptr-rsptr];
			++ptr;
		}

		// Getting all regions starting before the fragment end.
		ptr=ptr_copy;
		while (ptr!=rsptr) {
			--ptr;
			dist = curstart - *ptr;
			if (dist > maxrange) { break; }
			optr[dist < -maxrange ? -maxrange : dist] += wptr[ptr-rsptr];
			dist2 = curend - *ptr;
			if (dist2 < maxrange) { optr[dist2+1] -= wptr[ptr-rsptr]; }
		}
	}

	optr -= maxrange;
	for (int i=1; i<totallen; ++i) { optr[i]+=optr[i-1]; }
} catch (std::exception& e) {
	UNPROTECT(1);
	throw;
}

	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

