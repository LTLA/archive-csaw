#include "csaw.h"

SEXP get_cluster_stats (SEXP fcdex, SEXP pvaldex, SEXP tab, SEXP by, SEXP weight, SEXP fcthreshold) try {
	// Checking indices.
	if (!isInteger(fcdex) || !isInteger(pvaldex)) { throw std::runtime_error("table indices should be integer"); }
	if (LENGTH(pvaldex)!=1) { throw std::runtime_error("only one index should be supplied for the p-value column"); }
	const int pdex=asInteger(pvaldex);
	const int fcn=LENGTH(fcdex);
	const int* odptr=INTEGER(fcdex);

	// Setting up the p-value columns.
	if (!isNewList(tab)) { throw std::runtime_error("data values should be supplied as a list or dataframe"); }
    if (pdex < 0 || pdex >= LENGTH(tab)) { throw std::runtime_error("p-value index out of range"); }
	SEXP pval=VECTOR_ELT(tab, pdex);
	if (!isNumeric(pval)) { throw std::runtime_error("vector of p-values should be double precision"); }
	const double *pptr=REAL(pval);
	const int n=LENGTH(pval);

	// Setting up the log-FC columns.
	double** fcptrs=(double**)R_alloc(fcn, sizeof(double*));
	for (int i=0; i<fcn; ++i) { 
        const int& curfcdex=odptr[i];
        if (curfcdex < 0 || curfcdex >= LENGTH(tab)) { throw std::runtime_error("log-FC index out of range"); }
		SEXP logfc=VECTOR_ELT(tab, curfcdex);
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
	SEXP output=PROTECT(allocVector(VECSXP, 4));
	try {
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, total));
		SET_VECTOR_ELT(output, 1, allocMatrix(INTSXP, total, fcn*2));
		SET_VECTOR_ELT(output, 2, allocVector(REALSXP, total));
        SET_VECTOR_ELT(output, 3, allocVector(INTSXP, fcn==1 ? total : 0));
		if (total==0) {
			UNPROTECT(1);
			return output;
		}
        
		int* otptr=INTEGER(VECTOR_ELT(output, 0));
		int** ofptrs=(int**)R_alloc(fcn*2, sizeof(int*));
        if (fcn) {
            ofptrs[0]=INTEGER(VECTOR_ELT(output, 1));
            for (int i=1; i<fcn*2; ++i) { ofptrs[i]=ofptrs[i-1]+total; }
        }
		double* opptr=REAL(VECTOR_ELT(output, 2));
        int* odptr=INTEGER(VECTOR_ELT(output, 3));
	
        // Various temporary values.
		int i=0, j, k, x;
        double subweight, more_temp, remaining;
        int minx;
        bool has_up, has_down;
     
		while (i<n) {
			j=i+1;
			subweight=wptr[i];
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
			more_temp=0;
            remaining=wptr[sortvec[i]];
            minx=i;
			double& outp=(*opptr=pptr[sortvec[i]]/remaining); 

			for (x=i+1; x<j; ++x) {
 		    	const int& current=sortvec[x];	
				remaining+=wptr[current];
				more_temp=pptr[current]/remaining;
				if (more_temp < outp) { 
                    outp=more_temp; 
                    minx=x;
                }
			}
			outp*=subweight;

            /* If there's only one log-FC, we also determine which directions contribute to the combined p-value.
             * This is done by looking at the direction of the tests with p-values below that used as the combined p-value.
             * These tests must contribute because if any of them were non-significant, the combined p-value would increase.
             * Output codes are only up (1), only down (2) or mixed, i.e., both up and down (0).
             */
            if (fcn==1) {
                has_up=false;
                has_down=false;
                for (x=i; x<=minx; ++x) {
                    const double& curfc=fcptrs[0][sortvec[x]];
                    if (curfc > 0) { has_up=true; }
                    if (curfc < 0) { has_down=true; }
                    if (has_up & has_down) { break; }
                }
                int& current_dir=(odptr[0]=0);
                if (has_up & !has_down) { current_dir=1; }
                else if (has_down & !has_up) { current_dir=2; }
                ++odptr;
            }

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

/* Computes the total cluster weight in a reasonably fast manner. */

SEXP get_cluster_weight(SEXP ids, SEXP weight) try {
	if (!isInteger(ids)) { throw std::runtime_error("vector of cluster ids should be integer"); }
	if (!isNumeric(weight)) { throw std::runtime_error("vector of weights should be double precision"); }
	const double *wptr=REAL(weight);
	const int* iptr=INTEGER(ids);
    const int n=LENGTH(ids);
	if (n!=LENGTH(weight)) { throw std::runtime_error("vector lengths are not equal"); }

    int total=0;
	if (n > 0) {
		total=1;
		for (int i=1; i<n; ++i) {
			if (iptr[i] < iptr[i-1]) { throw std::runtime_error("vector of cluster ids should be sorted"); }
			else if (iptr[i]!=iptr[i-1]) { ++total; }
		}
	}

    SEXP output=PROTECT(allocVector(REALSXP, total));
	try {
        if (total) { 
            double* optr=REAL(output);
            (*optr)=wptr[0];
            for (int i=1; i<n; ++i) {
                if (iptr[i]!=iptr[i-1]) { 
                    ++optr; 
                    (*optr)=0;
                }
                (*optr)+=wptr[i];
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

