#include "csaw.h"

/* This function computes the mean and standard deviation of the count vectors
 * across the genome, at every delay distance. This accounts for the genomic
 * positions lost at the end of the chromosome (for the forward strand) or at
 * the start of the chromosome (on the reverse strand) when shifting occurs.
 *
 * The function operates by computing the mean and variance of the 'core' positions,
 * i.e., the positions that are present at the maximum delay distance. It then
 * computes the running mean and variance when the delay is wound back, using 
 * Welford's algorithm at every decreasing delay distance.
 * 
 * This approach is more numerically stable than the sum of products method to
 * compute the variance (and the standard variance calculation can't adjust easily
 * to changes in the mean). We also return the first delay distance which has a 
 * non-zero standard deviation, to avoid division by zero later.
 */

int fill (int n, std::deque<double>& mu, std::deque<double>& sd, const int* pos_ptr, const int* cnt_ptr,
		const int total_len, const int chr_len, const bool reverse) {
	mu.resize(n+1);
	sd.resize(n+1);
	if (n>chr_len-2) { n=chr_len-2; } // To ensure that 2 bases are available, so variance is estimable at max delay. 

	// Changing parameters according to forward/reverse reads.
	int step=1, start=0, end=total_len, threshold=chr_len-n;
	if (reverse) { 
		step=-1;
		start=total_len-1;
		end=-1;
		threshold=n+1;
	}

	// Find index of first read (going from end->start) with position <= threshold; adding 'step' to get endpoint for 'core' positions, i.e., from 'start' to 'temp_end'.
	int temp_end=start; 
	for (int i=end-step; i!=start-step; i-=step) {
		if ((threshold - pos_ptr[i])*step >= 0) { 
			temp_end=i+step;
			break;
		}
	}
	
	// Computing the mean over the 'core' positions.
	double& curmean=mu[n];
	for (int index=start; index!=temp_end; index+=step) { curmean += cnt_ptr[index]; }
	curmean/=chr_len-n;

	// Computing the variance for the 'core' positions. Using a jump table to speed things up with lots of zero's and one's (e.g. dedupped data).
	int first_pos_sd=-1;
	std::deque<bool> present;
	int present_counts=0;
	{
		double& curvar=sd[n];
		std::deque<double> jumptab(2, R_NaReal); // Setting it up to a reasonable extent (need at least 0 and 1).
		for (int index=start; index!=temp_end; index+=step) { 
			const int& curcount=cnt_ptr[index];
			if (int(jumptab.size()) <= curcount) { jumptab.resize(curcount+1, R_NaReal); } 
			if (ISNA(jumptab[curcount])) {
				jumptab[curcount]=curcount-curmean;
				jumptab[curcount]*=jumptab[curcount];
			}
			curvar += jumptab[curcount];
		}
		
		// Adding the squared product for all the positions with zeroes.
		const int num_zero_pos=chr_len-n-(temp_end-start)*step;
		curvar+=curmean*curmean*num_zero_pos; 

		// Checking whether it would theoretically have a zero variance, i.e., need at least two different integers in cnt_ptr (including implicit zeros). 
		present.resize(jumptab.size(), false);
		if (num_zero_pos) { 
			++present_counts;
			present[0]=true;
		}
		for (size_t i=1; i<jumptab.size(); ++i) { 
			if (!ISNA(jumptab[i])) { 
				++present_counts; 
				if (present_counts>=2) { 
					first_pos_sd=n; 
					break;
				}
				present[i]=true;
			} 
		}
	}
	
	// Computing the running mean and variance with decreasing delay distance, using Welford's algorithm as the next base position comes online.
	int current_pos=threshold,
		current_index=temp_end,
		current_count=0;
	double delta;
	for (int i=n-1; i>=0; --i) {
 		current_pos+=step;
		if (current_index!=end && pos_ptr[current_index]==current_pos) {
			current_count=cnt_ptr[current_index];
			current_index+=step;
		} else { current_count=0; }

		// Welford's algorithm [Technometrics, 4(3):419-420].
		delta=current_count-mu[i+1];
		mu[i]=mu[i+1]+delta/(chr_len-i);
		sd[i]=sd[i+1]+delta*(current_count-mu[i]);
		
		/* Checking if the variance is zero, and if so, whether the current base position gets a non-zero variance.
		 * Any addition will be sufficient, as 'counted' should be at least 1 by this point, if there is at least 1 read.
 	     * Note that this is 'first' in terms of decreasing delay distance; it'll be the largest delay distance for
		 * which the standard deviation is positive.
		 */
		if (first_pos_sd < 0) {
			if (current_count>=int(present.size())) { present.resize(current_count+1, false); }
			if (!present[current_count]) { 
				present[current_count]=true;
				++present_counts;
			}
			if (present_counts>2) { 
				throw std::runtime_error("first delay distance with positive s.d. should already be assigned"); 
			} else if (present_counts==2) {
				first_pos_sd=i; 
			}
		}
	}

	// Obtaining the standard deviation.
	for (int i=0; i<=n; ++i) {
		sd[i]/=chr_len-i-1;
		sd[i]=sqrt(sd[i]);
	}
	
	return first_pos_sd;
}

/* A function to calculate the correlations between reads on the same chromosome.
 * Cross-correlations can be calculated by supplying RLE value/lengths of forward reads 
 * as '1' and that of reverse reads as '2'. Autocorrelations can be calculated by 
 * supplying RLE value/lengths of all reads in both '1' and '2'.
 * 
 * It returns a vector of correlation coefficients from 0 (no shift) to the maximum shift.
 * We use the original expression for the correlation coefficient rather than using the
 * algebraically equivalent (but numerically unstable) sum of products method. There's
 * still a bit of numerical instability at the end where a subtraction step is involved,
 * but you can't have everything, I suppose.
 */

SEXP correlate_reads (SEXP pos1, SEXP num1, SEXP pos2, SEXP num2, SEXP max_dist, SEXP chrlen) try {
    // Prepping up the R input data.
    if (!isInteger(pos1)) { throw std::runtime_error("forward positions must be an integer vector"); }
    if (!isInteger(num1)) { throw std::runtime_error("forward counts must be an integer vector"); }
    if (!isInteger(pos2)) { throw std::runtime_error("reverse positions must be an integer vector"); }
    if (!isInteger(num2)) { throw std::runtime_error("reverse counts must be an integer vector"); }
    const int * fpptr=INTEGER(pos1);
    const int * rpptr=INTEGER(pos2);
    const int * fcptr=INTEGER(num1);
    const int * rcptr=INTEGER(num2);
    const int fLen=LENGTH(pos1), rLen=LENGTH(pos2);
    if (fLen!=LENGTH(num1) || rLen!=LENGTH(num2)) { 
       	throw std::runtime_error("lengths of position vectors do not correspond to frequency vectors"); }

	// Checking other scalars.
    if (!isInteger(max_dist) || LENGTH(max_dist)!=1) { throw std::runtime_error("maximum distance should be an integer scalar"); }	
    const int mdist=asInteger(max_dist);
	if (mdist <= 0) { throw std::runtime_error("maximum distance should be a positive integer"); }
	if (!isInteger(chrlen) || LENGTH(chrlen)!=1) { throw std::runtime_error("length of chromosome must be an integer scalar");} 
	const int clen=asInteger(chrlen);

	// Computing the mean and variance.
	std::deque<double> fmean, rmean, fsd, rsd;
	const int ffirst=fill(mdist, fmean, fsd, fpptr, fcptr, fLen, clen, false);
	const int rfirst=fill(mdist, rmean, rsd, rpptr, rcptr, rLen, clen, true);

    // Setting up to go through the forwards.
    SEXP sumOut=PROTECT(allocVector(REALSXP, mdist+1));
	try {
    	double* sumptr=REAL(sumOut);
    	for (int i=0; i<=mdist; ++i) { sumptr[i]=0; }
		std::deque<double> sumfdiff(mdist+1), sumrdiff(mdist+1);
		std::deque<int> nonempty(mdist+1);

    	// Looking at every forward/reverse combination, where the reverse read must be downstream of the forward read.
		double rdiff, fdiff;
    	int f2rdex=0, f2rdex_copy;
    	for (int fdex=0; fdex<fLen; ++fdex) {
        	const int& fpos=fpptr[fdex];
        	const int& fcounts=fcptr[fdex];

        	// Going through the reverse strand.
        	f2rdex_copy=f2rdex;
        	while (f2rdex_copy < rLen) {
            	const int& rpos=rpptr[f2rdex_copy];
            	const int diff=rpos-fpos;
            	if (diff < 0) {
                	++f2rdex;
            	} else if (diff > mdist) { 
					break; 
				} else {
					fdiff=fcounts-fmean[diff];
					rdiff=rcptr[f2rdex_copy]-rmean[diff];
                	sumptr[diff]+=rdiff*fdiff;
					sumfdiff[diff]+=fdiff;
					sumrdiff[diff]+=rdiff;
					++nonempty[diff];
            	}
            	++f2rdex_copy; 
        	}
    	}

		/* Computing the actual correlation from the various collected statistics. The correlation is defined as:
		 *     sum of PRODUCT=(x-_X)(y-_Y) for all x and y
		 *   = SUM[PRODUCT] for x>0,y>0 + SUM[PRODUCT] for x==0,y>0 + SUM[PRODUCT] for x>0,y==0 + SUM[PRODUCT] for x==0,y==0
		 *
 		 * The first expression is that in sum_ptr. The second is equal to 
		 *     -_X*( SUM[y-_Y] for x==0,y>0 )
 		 *   = -_X*( SUM[y-_Y] for all x, y - SUM[y-_Y] for x>0,y>0 - SUM[y-_Y] for x==0,y==0 - SUM[y-_Y] for x>0,y==0)
 		 *   = -_X*( 0 - SUM[y-_Y] for x>0,y>0 - SUM[y-_Y] for x==0,y==0 - SUM[y-_Y] for x>0,y==0)
 		 *   = -_X*( 0 - SUM[y-_Y] for x>0,y>0 + _Y*( # of times y==0 ) )
		 *
 		 * where the first term inside the brackets is simply that of sumrdiff. The third expression, by similarity, is that of 
		 *      SUM[PRODUCT] for x>0,y==0
 		 *    = -_Y*( - SUM[x-_X] for x>0,y>0 + _X*( # of times x==0 ) )
		 *
 		 * where the first term inside the brackets is simply that of sumfdiff. And the last expression is simply 
		 *      _X*_Y*( # of times x==0 and y==0 ).
 		 *
 		 * When assembled, we get:
  	   	 *      sum_ptr[i] + _X*sumrdiff[i] + _Y*sumfdiff[i] - _X*_Y*( # of times x==0 + # of times y==0 - # of times both==0)
		 *
 		 * The expression inside the brackets is simply the number of pairs where either is equal to zero; which, in turn, is equal 
		 * to the total number of pairs (after shifting by the delay) minus the number of pairs where both are non-zero.
 		 *
 		 * Subtraction of the final term could be numerically unstable. This can't really be helped if speed is to be maintained.
		 * In practice, it'll probably be okay because rmean and fmean should be fairly small (a lot more base positions than reads).
 		 */
		for (int i=0; i<=mdist && i<=clen-2; ++i) {
			// Providing some protection for delay distances where the s.d. is zero.
			if (i>ffirst || i>rfirst) {
				sumptr[i]=0;
			} else {
				sumptr[i]-=rmean[i]*fmean[i]*(clen-i-nonempty[i]); // A bit of numerical instability with this step.
				sumptr[i]+=rmean[i]*sumfdiff[i]+fmean[i]*sumrdiff[i];
				sumptr[i]/=rsd[i]*fsd[i];
				sumptr[i]/=clen-i-1; 
			}
		}
  	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}	  
    
	UNPROTECT(1);
    return sumOut;
} catch (std::exception& e) {
	return mkString(e.what()); 
}
