#include "csaw.h"

extern "C" {


SEXP R_get_cluster_stats (SEXP fcdex, SEXP cpmdex, SEXP pvaldex, SEXP tab, SEXP by, SEXP weight) try {
	// Checking indices.
	if (!IS_INTEGER(fcdex) || !IS_INTEGER(cpmdex) || !IS_INTEGER(pvaldex)) { throw std::runtime_error("table indices should be integer"); }
	if (LENGTH(cpmdex)!=1 || LENGTH(pvaldex)!=1) { throw std::runtime_error("only one index should be supplied for log-CPM and p-value columns"); }
	const int pdex=INTEGER_VALUE(pvaldex), cdex=INTEGER_VALUE(cpmdex);
	const int fcn=LENGTH(fcdex);
	if (!LENGTH(fcdex)) { throw std::runtime_error("at least one index should be supplied for log-FC columns"); }
	const int* fcdptr=INTEGER_POINTER(fcdex);

	// Setting up the columns.
	if (!isNewList(tab)) { throw std::runtime_error("data values should be supplied as a list or dataframe"); }
	SEXP cpm=VECTOR_ELT(tab, cdex), pval=VECTOR_ELT(tab, pdex);
	if (!IS_NUMERIC(cpm)) { throw std::runtime_error("vector of log-CPMs should be double precision"); }
	if (!IS_NUMERIC(pval)) { throw std::runtime_error("vector of p-values should be double precision"); }
	const double *cptr=NUMERIC_POINTER(cpm), *pptr=NUMERIC_POINTER(pval);
	const int n=LENGTH(pval);
	if (n!=LENGTH(cpm)) { throw std::runtime_error("vector lengths are not equal"); }

	// Setting up the columns (II), for log-FCs.
	double** fcptrs=(double**)R_alloc(fcn, sizeof(double*));
	for (int i=0; i<fcn; ++i) { 
		SEXP fc=VECTOR_ELT(tab, fcdptr[i]);
		if (!IS_NUMERIC(fc)) { throw std::runtime_error("vector of log-FCs should be double precision"); }
	    if (n!=LENGTH(fc)) { throw std::runtime_error("vector lengths are not equal"); }
		fcptrs[i]=NUMERIC_POINTER(fc);
	}

	// Setting up the remaining inputs. 
	if (!IS_INTEGER(by)) { throw std::runtime_error("vector of cluster ids should be integer"); }
	if (!IS_NUMERIC(weight)) { throw std::runtime_error("vector of weights should be double precision"); }
	const double *wptr=NUMERIC_POINTER(weight);
	const int* bptr=INTEGER_POINTER(by);
	if (!n) { throw std::runtime_error("nothing to cluster"); }
	if (n!=LENGTH(by) || n!=LENGTH(weight)) { throw std::runtime_error("vector lengths are not equal"); }

	// Checking that the 'by' is sorted, counting the number of elements and setting up a vector of [0, n).
	int total=1;
	int* sortvec=(int*)R_alloc(n, sizeof(int));
	sortvec[0]=0;
	for (int i=1; i<n; ++i) { 
		if (bptr[i] < bptr[i-1]) { throw std::runtime_error("vector of cluster ids should be sorted"); }
		else if (bptr[i]!=bptr[i-1]) { ++total; }
		sortvec[i]=i;
	}
	sort_row_index<double> pcomp(pptr);

	// Pulling out results.
	SEXP output=PROTECT(NEW_LIST(3));
	try {
		SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, total, fcn));
		double** ofptrs=(double**)R_alloc(fcn, sizeof(double*));
		ofptrs[0]=NUMERIC_POINTER(VECTOR_ELT(output, 0));
		for (int i=1; i<fcn; ++i) { ofptrs[i]=ofptrs[i-1]+total; }
		SET_VECTOR_ELT(output, 1, NEW_NUMERIC(total));
		double* ocptr=NUMERIC_POINTER(VECTOR_ELT(output, 1));
		SET_VECTOR_ELT(output, 2, NEW_NUMERIC(total));
		double* opptr=NUMERIC_POINTER(VECTOR_ELT(output, 2));
	
		int i=0;
		while (i<n) {
			int j=i+1;
			double subweight=wptr[i];
			while (j < n && bptr[i]==bptr[j]) { 
				subweight+=wptr[j];
				++j; 
			}

			// Computing the weighted mean of log-FC(s) and log-CPM values.
			double& outcpm=(*ocptr=0);
	        for (int x=i; x<j; ++x) { outcpm+=cptr[x]*wptr[x]; }
			for (int k=0; k<fcn; ++k) { 
				double& outfc=(ofptrs[k][0]=0);
				for (int x=i; x<j; ++x) {  outfc+=fcptrs[k][x]*wptr[x]; }
				outfc/=subweight;
				++(ofptrs[k]);
			}
			outcpm/=subweight;

			/* Computing the weighted Simes value. The weights are implemented as frequency 
 			 * weights, e.g., if you had 2 tests with a weight of 10 to 1, you'd consider the
 			 * one with the higher weight 10 more times to try to reject the global null (i.e.,
 			 * expanding it in-place in the sorted vector of p-values).
 			 */
			std::sort(sortvec+i, sortvec+j, pcomp);
			double more_temp=0, remaining=wptr[sortvec[i]];
			double& outp=(*opptr=pptr[sortvec[i]]/remaining); 
			for (int x=i+1; x<j; ++x) {
 		    	const int& current=sortvec[x];	
				remaining+=wptr[current];
				more_temp=pptr[current]/remaining;
				if (more_temp < outp) { outp=more_temp; }
			}
			outp*=subweight;

			// Setting it up for the next round.
			++ocptr;
			++opptr;
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

