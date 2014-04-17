#include "csaw.h"

extern "C" {

SEXP R_best_in_cluster(SEXP pval, SEXP by, SEXP weight) try {
	if (!IS_NUMERIC(pval)) { throw std::runtime_error("vector of p-values should be double precision"); }
	const double *pptr=NUMERIC_POINTER(pval);
	const int n=LENGTH(pval);
	if (!IS_INTEGER(by)) { throw std::runtime_error("vector of cluster ids should be integer"); }
	if (!IS_NUMERIC(weight)) { throw std::runtime_error("vector of weights should be double precision"); }
	const double *wptr=NUMERIC_POINTER(weight);
	const int* bptr=INTEGER_POINTER(by);
	if (!n) { throw std::runtime_error("nothing to cluster"); }
	if (n!=LENGTH(by) || n!=LENGTH(weight)) { throw std::runtime_error("vector lengths are not equal"); }

	// Checking that the 'by' is sorted, counting the number of elements.
	int total=1;
	for (int i=1; i<n; ++i) { 
		if (bptr[i] < bptr[i-1]) { throw std::runtime_error("vector of cluster ids should be sorted"); }
		else if (bptr[i]!=bptr[i-1]) { ++total; }
	}

	// Pulling out results.
	SEXP output=PROTECT(NEW_LIST(2));
	try {
		SET_VECTOR_ELT(output, 0, NEW_NUMERIC(total));
		double* opptr=NUMERIC_POINTER(VECTOR_ELT(output, 0));
		SET_VECTOR_ELT(output, 1, NEW_INTEGER(total));
		int* oiptr=INTEGER_POINTER(VECTOR_ELT(output, 1));
	
		int i=0;
		while (i<n) {
			int j=i+1;
			double subweight=wptr[i];
			while (j < n && bptr[i]==bptr[j]) { 
				subweight+=wptr[j];
				++j; 
			}

			/* Computing the Holm p-value for the best window (basically Bonferroni, if we're taking the minimum).
 			 * Weights are considered as relative frequency weights i.e. the total number of tests is rescaled
 			 * relative to the weight of the current test (so, [10,1] weights would consider there to be 1.1 tests
 			 * for the first one and 11 tests for the second one).
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

}

