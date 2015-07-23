#include "csaw.h"

SEXP get_cluster_stats (SEXP fcdex, SEXP pvaldex, SEXP tab, SEXP by, SEXP weight, SEXP fcthreshold) try {
	// Checking indices.
	if (!isInteger(fcdex) || !isInteger(pvaldex)) { throw std::runtime_error("table indices should be integer"); }
	if (LENGTH(pvaldex)!=1) { throw std::runtime_error("only one index should be supplied for the p-value column"); }
	const int pdex=asInteger(pvaldex);
	const int fcn=LENGTH(fcdex);
	if (!fcn) { throw std::runtime_error("at least one index should be supplied for the log-FC columns"); }
	const int* odptr=INTEGER(fcdex);

	// Setting up the p-value columns.
	if (!isNewList(tab)) { throw std::runtime_error("data values should be supplied as a list or dataframe"); }
	SEXP pval=VECTOR_ELT(tab, pdex);
	if (!isNumeric(pval)) { throw std::runtime_error("vector of p-values should be double precision"); }
	const double *pptr=REAL(pval);
	const int n=LENGTH(pval);

	// Setting up the log-FC columns.
	double** fcptrs=(double**)R_alloc(fcn, sizeof(double*));
	for (int i=0; i<fcn; ++i) { 
		SEXP logfc=VECTOR_ELT(tab, odptr[i]);
		if (!isNumeric(logfc)) { throw std::runtime_error("vector of logfc statistics should be double precision"); }
		if (n!=LENGTH(logfc)) { throw std::runtime_error("vector lengths are not equal"); }
		fcptrs[i]=REAL(logfc);
	}
	if (!isReal(fcthreshold) || LENGTH(fcthreshold)!=1) { throw std::runtime_error("log-fold change threshold should be a numeric scalar"); }
	const double fcthresh=asReal(fcthreshold);

	// Setting up the remaining inputs. 
	if (!isInteger(by)) { throw std::runtime_error("vector of cluster ids should be integer"); }
	if (!isNumeric(weight)) { throw std::runtime_error("vector of weights should be double precision"); }
	const double *wptr=REAL(weight);
	const int* bptr=INTEGER(by);
	if (n!=LENGTH(by) || n!=LENGTH(weight)) { throw std::runtime_error("vector lengths are not equal"); }

	// Checking that the 'by' is sorted, counting the number of elements and setting up a vector of [0, n).
	int total=0;
	int* sortvec=NULL;
	sort_row_index<double> pcomp(pptr);
	if (n > 0) {
		total=1;
		sortvec=(int*)R_alloc(n, sizeof(int));
		sortvec[0]=0;
		for (int i=1; i<n; ++i) {
			if (bptr[i] < bptr[i-1]) { throw std::runtime_error("vector of cluster ids should be sorted"); }
			else if (bptr[i]!=bptr[i-1]) { ++total; }
			sortvec[i]=i;
		}
	}

	// Pulling out results.
	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, total));
		SET_VECTOR_ELT(output, 1, allocMatrix(INTSXP, total, fcn*2));
		SET_VECTOR_ELT(output, 2, allocVector(REALSXP, total));
		if (total==0) {
			UNPROTECT(1);
			return output;
		}

		int* otptr=INTEGER(VECTOR_ELT(output, 0));
		int** ofptrs=(int**)R_alloc(fcn, sizeof(int*));
		ofptrs[0]=INTEGER(VECTOR_ELT(output, 1));
		for (int i=1; i<fcn*2; ++i) { ofptrs[i]=ofptrs[i-1]+total; }
		double* opptr=REAL(VECTOR_ELT(output, 2));
	
		int i=0, k, x;
		while (i<n) {
			int j=i+1;
			double subweight=wptr[i];
			while (j < n && bptr[i]==bptr[j]) { 
				subweight+=wptr[j];
				++j; 
			}

			// Computing the total number of windows, and that up/down for each fold change.
			*(otptr)=j-i;
			++otptr;
			for (k=0; k<fcn; ++k) { 
				int& allup=(ofptrs[k*2][0]=0);
				int& alldown=(ofptrs[k*2+1][0]=0);
				for (x=i; x<j; ++x) {  
					if (fcptrs[k][x] > fcthresh) { ++allup; } 
					else if (fcptrs[k][x] < -fcthresh) { ++alldown; }
				}
				++(ofptrs[k*2]);
				++(ofptrs[k*2+1]);
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
