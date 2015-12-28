#include "csaw.h"

SEXP best_in_cluster(SEXP pval, SEXP by, SEXP weight) try {
	if (!isNumeric(pval)) { throw std::runtime_error("vector of p-values should be double precision"); }
	const double *pptr=REAL(pval);
	const int n=LENGTH(pval);
	if (!isInteger(by)) { throw std::runtime_error("vector of cluster ids should be integer"); }
	if (!isNumeric(weight)) { throw std::runtime_error("vector of weights should be double precision"); }
	const double *wptr=REAL(weight);
	const int* bptr=INTEGER(by);
	if (n!=LENGTH(by) || n!=LENGTH(weight)) { throw std::runtime_error("vector lengths are not equal"); }

	// Checking that the 'by' is sorted, counting the number of elements.
	int total=0;
	if (n > 0) {
		total=1;
		for (int i=1; i<n; ++i) {
			if (bptr[i] < bptr[i-1]) { throw std::runtime_error("vector of cluster ids should be sorted"); }
			else if (bptr[i]!=bptr[i-1]) { ++total; }
		}
	}

	// Pulling out results.
	SEXP output=PROTECT(allocVector(VECSXP, 2));
	try {
		SET_VECTOR_ELT(output, 0, allocVector(REALSXP, total));
		SET_VECTOR_ELT(output, 1, allocVector(INTSXP, total));
		if (total==0) {
			UNPROTECT(1);
			return output;
		}
		double* opptr=REAL(VECTOR_ELT(output, 0));
		int* oiptr=INTEGER(VECTOR_ELT(output, 1));
	
		int i=0;
		while (i<n) {
			int j=i+1;
			double subweight=wptr[i];
			while (j < n && bptr[i]==bptr[j]) { 
				subweight+=wptr[j];
				++j; 
			}

			/* Computing the Holm p-value for the best window (basically Bonferroni, if we're taking the minimum).
			 * Weights are defined according to the weighted Bonferroni (see http://arxiv.org/abs/math.ST/0604172,
			 * though some mental arithmetic is needed). These can also be treated as relative frequency weights,
			 * i.e. the total number of tests is rescaled relative to the weight of the current test (so, [10,1] 
			 * weights would consider there to be 1.1 tests for the first one and 11 tests for the second one).
			 */
			int& outi=(*oiptr=i);
			double& outp=(*opptr=pptr[outi]/wptr[outi]);
			double tempp=0;
			for (int x=i+1; x<j; ++x) {
				tempp=pptr[x]/wptr[x];
				if (tempp < outp) { 
					outi=x;
					outp=tempp;
				}
			}
			outp*=subweight;
 		    if (outp > 1) { outp=1; }	
			++outi; // For 1-based indexing in R.

			// Setting it up for the next round.
			++opptr;
			++oiptr;
			i=j;
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
