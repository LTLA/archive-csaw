#include "csaw.h"

SEXP get_cluster_stats (SEXP otherdex, SEXP pvaldex, SEXP tab, SEXP by, SEXP weight) try {
	// Checking indices.
	if (!isInteger(otherdex) || !isInteger(pvaldex)) { throw std::runtime_error("table indices should be integer"); }
	if (LENGTH(pvaldex)!=1) { throw std::runtime_error("only one index should be supplied for log-CPM and p-value columns"); }
	const int pdex=asInteger(pvaldex);
	const int ocn=LENGTH(otherdex);
	if (!ocn) { throw std::runtime_error("at least one index should be supplied for the other columns"); }
	const int* odptr=INTEGER(otherdex);

	// Setting up the columns.
	if (!isNewList(tab)) { throw std::runtime_error("data values should be supplied as a list or dataframe"); }
	SEXP pval=VECTOR_ELT(tab, pdex);
	if (!isNumeric(pval)) { throw std::runtime_error("vector of p-values should be double precision"); }
	const double *pptr=REAL(pval);
	const int n=LENGTH(pval);
	if (n==0) { throw std::runtime_error("no elements supplied to compute cluster statistics"); }

	// Setting up the columns (II), for log-FCs.
	double** optrs=(double**)R_alloc(ocn, sizeof(double*));
	for (int i=0; i<ocn; ++i) { 
		SEXP other=VECTOR_ELT(tab, odptr[i]);
		if (!isNumeric(other)) { throw std::runtime_error("vector of other statistics should be double precision"); }
	    if (n!=LENGTH(other)) { throw std::runtime_error("vector lengths are not equal"); }
		optrs[i]=REAL(other);
	}

	// Setting up the remaining inputs. 
	if (!isInteger(by)) { throw std::runtime_error("vector of cluster ids should be integer"); }
	if (!isNumeric(weight)) { throw std::runtime_error("vector of weights should be double precision"); }
	const double *wptr=REAL(weight);
	const int* bptr=INTEGER(by);
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
	SEXP output=PROTECT(allocVector(VECSXP, 2));
	try {
		SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, total, ocn));
		double** ooptrs=(double**)R_alloc(ocn, sizeof(double*));
		ooptrs[0]=REAL(VECTOR_ELT(output, 0));
		for (int i=1; i<ocn; ++i) { ooptrs[i]=ooptrs[i-1]+total; }
		SET_VECTOR_ELT(output, 1, allocVector(REALSXP, total));
		double* opptr=REAL(VECTOR_ELT(output, 1));
	
		int i=0, k, x;
		while (i<n) {
			int j=i+1;
			double subweight=wptr[i];
			while (j < n && bptr[i]==bptr[j]) { 
				subweight+=wptr[j];
				++j; 
			}

			// Computing the weighted mean of log-FC(s) and log-CPM values.
			for (k=0; k<ocn; ++k) { 
				double& out_other=(ooptrs[k][0]=0);
				for (x=i; x<j; ++x) {  out_other+=optrs[k][x]*wptr[x]; }
				out_other/=subweight;
				++(ooptrs[k]);
			}

			/* Computing the weighted Simes value. The weights are implemented as frequency 
 			 * weights, e.g., if you had 2 tests with a weight of 10 to 1, you'd consider the
 			 * one with the higher weight 10 more times to try to reject the global null (i.e.,
 			 * expanding it in-place in the sorted vector of p-values).
 			 */
			std::sort(sortvec+i, sortvec+j, pcomp);
			double more_temp=0, remaining=wptr[sortvec[i]];
			double& outp=(*opptr=pptr[sortvec[i]]/remaining); 
			for (x=i+1; x<j; ++x) {
 		    	const int& current=sortvec[x];	
				remaining+=wptr[current];
				more_temp=pptr[current]/remaining;
				if (more_temp < outp) { outp=more_temp; }
			}
			outp*=subweight;

			// Setting it up for the next round.
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
