#include "csaw.h"

#ifdef DEBUG
#include <map>
#endif

/* This function scans through the track and pulls out local maxima. */

SEXP get_profile(SEXP starts, SEXP ends, SEXP regstarts, SEXP regends, 
		SEXP order_s, SEXP order_e, SEXP rank_e, SEXP weights, SEXP range) try {

	if (!isInteger(starts) || !isInteger(ends)) { throw std::runtime_error("fragment start/end positions should be integer"); }
	const int nfrags=LENGTH(starts);
	if (LENGTH(ends)!=nfrags) { throw std::runtime_error("fragment start/end vectors should have same length"); }
	if (!isInteger(regstarts) || !isInteger(regends)) { throw std::runtime_error("region start/end positions should be integer"); }
	const int nregs=LENGTH(regstarts);
	if (LENGTH(regends)!=nregs) { throw std::runtime_error("region start/end vectors should have the same length"); }
	if (nregs==0) { throw std::runtime_error("no regions supplied"); }

	if (!isInteger(order_s) || !isInteger(order_e) || !isInteger(rank_e)) { throw std::runtime_error("ordering vectors should be integer"); }
	if (LENGTH(order_s)!=nregs || LENGTH(order_e)!=nregs || LENGTH(rank_e)!=nregs) { 
		throw std::runtime_error("ordering vectors should have length equal to number of regions"); }
	if (!isReal(weights))  { throw std::runtime_error("weight vector should be double-precision"); }
	if (LENGTH(weights)!=nregs) { throw std::runtime_error("weight vector should have length equal to number of regions"); }

	// Setting up constructs.
	if (!isInteger(range) || LENGTH(range)!=1) { throw std::runtime_error("range distance should be an integer scalar"); }
	const int maxrange=asInteger(range);
	const int* fsptr=INTEGER(starts),
		  *feptr=INTEGER(ends),
		  *rsptr=INTEGER(regstarts),
		  *reptr=INTEGER(regends),
		  *osptr=INTEGER(order_s),
		  *oeptr=INTEGER(order_e),
		  *erptr=INTEGER(rank_e);
	const double* wptr=REAL(weights);

	// Setting up output.
	SEXP output=PROTECT(allocVector(REALSXP, maxrange));
try{ 
	double* optr=REAL(output);
	for (int i=0; i<maxrange; ++i) { optr[i]=0; }

	// Running through the reads.
	const int* ptr=NULL, *ptr_copy=NULL, *terminal_s=rsptr+nregs;
	int dist, counter, end_rank;
	for (int frag=0; frag<nfrags; ++frag) {

		// Getting all regions starting after the fragment.
		Rprintf("Aggregating start:\n");
		const int& curend=feptr[frag];
		ptr_copy=ptr=std::upper_bound(rsptr, terminal_s, curend);
		while (ptr!=terminal_s) { 
			dist = *ptr - curend - 1;
			if (dist >= maxrange) { break; }
			optr[dist] += 1/wptr[osptr[ptr-rsptr]];
			++ptr;
		}

		// Getting the first region ending before the fragment.
		Rprintf("Jumping to end:\n");
		counter=ptr_copy-rsptr;
		if (counter==nregs) { --counter; }		
		end_rank=erptr[osptr[counter]];
		const int& curstart=fsptr[frag];
		while (curstart <= reptr[end_rank]) {
			--end_rank;
			if (end_rank < 0) { break; }				
		}

		// Adding all regions ending before the fragment.
		Rprintf("Aggregating ends:\n");
		while (end_rank >= 0) {
			dist = curstart - reptr[end_rank] - 1;
			if (dist >= maxrange) { break; }
			optr[dist] += 1/wptr[oeptr[end_rank]];
			--end_rank;
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

