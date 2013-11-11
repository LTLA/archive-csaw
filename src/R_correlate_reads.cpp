#include "csaw.h"

/* Stable running calculation of the variance and mean, using Welford's algorithm for both values
 * at decreasing delay distances (i.e. calculating values for the common regions at the maximum
 * shift, then iterating across the delay distances and computing running values). This approach
 * is more numerically stable than the sum of products method (see below).
 *
 * We also return the first delay distance which has a positive standard deviation. This is easier
 * than doing it outside, where we have to tangle with numerical imprecision.
 */

int fill (int n, std::deque<double>& mu, std::deque<double>& sd, const int* pos_ptr, const int* cnt_ptr,
		const int total_len, const int chr_len, const bool reverse) {
	/* Generating vectors. Also avoiding issues where n>=chr_len, by resetting the maximum delay
 	 * shift so that there's at least two elements. Otherwise, the variance is undefined.
 	 */
	mu.resize(n+1);
	sd.resize(n+1);
	if (n>=chr_len) { n=chr_len-2; }
	int first_pos_sd=-1;

	int step=1, start=0, end=total_len, threshold=chr_len-n;
	if (reverse) { 
		step=-1;
		start=total_len-1;
		end=-1;
		threshold=n+1;
	}
	/* Treading backwards, and identifying the first index where the position
 	 * is past the threshold (in terms of increasing step, so we've got to +step
 	 * after we reach the point of equality).
 	 */
	int temp_end=start;
	for (int i=end-step; i!=start-step; i-=step) {
		if ((threshold - pos_ptr[i])*step >= 0) { 
			temp_end=i+step;
			break;
		}
	}

	/* Computing the mean and variance for the longest common subsection (i.e. the largest delay).
 	 * We use a jump table for the latter, just to speed things up when we're dealing with a large
 	 * number of zero's and one's.
 	 */
	double& curmean=mu[n];
	for (int index=start; index!=temp_end; index+=step) { curmean += cnt_ptr[index]; }
	curmean/=chr_len-n;
	double& curvar=sd[n];
	std::deque<bool> present;
	{
		std::deque<double> tmps(1, -1);
		for (int index=start; index!=temp_end; index+=step) { 
			const int& curcount=cnt_ptr[index];
			if (tmps.size() <= curcount) { tmps.resize(curcount+1, -1); } 
			if (tmps[curcount]<0) {
				tmps[curcount]=curcount-curmean;
				tmps[curcount]*=tmps[curcount];
			}
			curvar += tmps[curcount];
		}
		const int leftover=chr_len-n-(temp_end-start)*step;
		curvar+=curmean*curmean*leftover;

		/* Checking whether it would theoretically have a zero-sd, using integers. If there's 
 		 * at least two unique positive integers, or there's one unique positive integer and at 
 		 * least one zero value, then the standard deviation must be positive.
 		 */
		present.resize(tmps.size(), false);
		int counted=0;
		if (leftover) { 
			++counted;
			present[0]=true;
		}
		for (int i=1; i<tmps.size(); ++i) { 
			if (tmps[i]>0) { 
				++counted; 
				if (counted>=2) { 
					first_pos_sd=n; 
					break;
				}
				present[i]=true;
			} 
		}
	}
	
	/* Now, computing running mean and variance values for all the delays. We use Welford's algorithm
	 * to provide some numerical stability. We take the square root to get the standard deviation.
	 */
	int running = chr_len-n;
	int current_pos=threshold;
	int current_count=0;
	double delta;
	for (int i=n-1; i>=0; --i) {
		running+=1;
		current_pos+=step;
		if (pos_ptr[temp_end]==current_pos) {
			current_count=cnt_ptr[temp_end];
			temp_end+=step;
		} else { current_count=0; }

		delta=current_count-mu[i+1];
		mu[i]=mu[i+1]+delta/running;
		sd[i]=sd[i+1]+delta*(current_count-mu[i]);
		sd[i+1]/=running-2;
		sd[i+1]=sqrt(sd[i+1]);

		/* Also checking whether this would break the zero-sd condition. Any addition
		 * will be sufficient, as 'counted' should be at least 1 (i.e. at least two
		 * non-negative values should be in there, somewhere).
		 */
		if (first_pos_sd < 0 && (current_count>=present.size() || !present[current_count])) { 
			first_pos_sd=i; }
	}
	sd[0]/=chr_len-1;
	sd[0]=sqrt(sd[0]);

	return first_pos_sd;
}

/* A function to calculate the correlations between reads on the same chromosome.
 * Cross correlations can be calculated by supplying RLE value/lengths of forward reads 
 * as '1' and that of reverse reads as '2'. Autocorrelations can be calculated by 
 * supplying RLE value/lengths of all reads in both '1' and '2'.
 * 
 * It returns a vector of correlation coefficients from 0 (no shift) to the maximum shift.
 * We use the original expression for the correlation coefficient rather than using the
 * algebraically equivalent (but numerically unstable) sum of products method. There's
 * still a bit of numerical instability at the end where a subtraction step is involved,
 * but you can't have everything, I suppose.
 */

extern "C" {

SEXP R_correlate_reads (SEXP pos1, SEXP num1, SEXP pos2, SEXP num2, SEXP max_dist, SEXP chrlen) try {
    // Prepping up the R input data.
    if (!IS_INTEGER(pos1)) { throw std::runtime_error("forward positions must be an integer vector"); }
    if (!IS_INTEGER(num1)) { throw std::runtime_error("forward counts must be an integer vector"); }
    if (!IS_INTEGER(pos2)) { throw std::runtime_error("reverse positions must be an integer vector"); }
    if (!IS_INTEGER(num2)) { throw std::runtime_error("reverse counts must be an integer vector"); }
    const int * fpptr=INTEGER_POINTER(pos1);
    const int * rpptr=INTEGER_POINTER(pos2);
    const int * fcptr=INTEGER_POINTER(num1);
    const int * rcptr=INTEGER_POINTER(num2);
    const int fLen=LENGTH(pos1), rLen=LENGTH(pos2);
    if (fLen!=LENGTH(num1) || rLen!=LENGTH(num2)) { 
       	throw std::runtime_error("lengths of position vectors do not correspond to frequency vectors"); }

	// Checking other scalars.
    if (!IS_INTEGER(max_dist) || LENGTH(max_dist)!=1) { throw std::runtime_error("maximum distance should be an integer scalar"); }	
    const int mdist=INTEGER_VALUE(max_dist);
	if (mdist <= 0) { throw std::runtime_error("maximum distance should be a positive integer"); }
	if (!IS_INTEGER(chrlen) || LENGTH(chrlen)!=1) { throw std::runtime_error("length of chromosome must be an integer scalar");} 
	const int clen=INTEGER_VALUE(chrlen);

	// Computing the mean and variance.
	std::deque<double> fmean, rmean, fsd, rsd;
	const int ffirst=fill(mdist, fmean, fsd, fpptr, fcptr, fLen, clen, false);
	const int rfirst=fill(mdist, rmean, rsd, rpptr, rcptr, rLen, clen, true);

    // Setting up to go through the forwards.
    SEXP sumOut=PROTECT(allocVector(REALSXP, mdist+1));
	try {
    	double* sumptr=NUMERIC_POINTER(sumOut);
    	for (int i=0; i<=mdist; ++i) { sumptr[i]=0; }
		std::deque<double> sumfdiff(mdist+1), sumrdiff(mdist+1);
		std::deque<int> nonempty(mdist+1);

    	/* Now iterating through the valid forward/reverse combinations
    	 * (where the reverse read must be downstream of the forward read).
    	 * This is more efficient as we only have to consider the combinations present.
    	 */
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

		/* Computing the actual correlation from the various collected statistics. The correlation is 
 		 * defined as the sum of PRODUCT=(x-_X)(y-_Y) for all x and y. This can be written as 
 		 * SUM[PRODUCT] for x>0,y>0 + SUM[PRODUCT] for x==0,y>0 + SUM[PRODUCT] for x>0,y==0 + SUM[PRODUCT] for x==0,y==0.
 		 * The first expression is that in sum_ptr. The second is equal to -_X*( SUM[y-_Y] for x==0,y>0 )
 		 *    = -_X*( -SUM[y-_Y] for x>0,y>0 - SUM[y-_Y] for x==0,y==0 - SUM[y-_Y] for x==0,y>0 )
 		 *    = -_X*( -SUM[y-_Y] for x>0,y>0 - _Y*( # of times x==0 ) )
 		 *  	where the first term inside the brackets is simply that of sumrdiff.
 		 * The third expression, by similarity, is that of SUM[PRODUCT] for x>0,y==0
 		 *    = -_Y*( -SUM[y-_Y] for x>0,y>0 - _X*( # of times y==0 ) )
 		 *      where the first term inside the brackets is simply that of sumfdiff.
 		 * And the last expression is simply _X*_Y*( # of times x==0 and y==0 ).
 		 *
 		 * When assembled, we get something like sum_ptr[i] + _X*sumrdiff[i] + _Y*sumfdiff[i]
 		 *      - _X*_Y*( # of times x==0 + # of times y==0 - # of times both==0)
 		 * The expression inside the brackets is simply the number of pairs
 		 * where either is equal to zero; which, in turn, is equal to the total
 		 * number of pairs without the number of pairs where both are non-zero.
 		 *
 		 * The final subtraction is a bit numerically unstable, but if it done before the addition
 		 * of '_X*sumrdiff[i] + _Y*sumfdiff[i]', it should be equivalent to straight summation
 		 * where all x>0,y>0 pairs are in one stretch; all x==0,y==0 pairs in the next stretch; and
 		 * all other pairs in the final stretch. It's not quite the worst case where all positive
 		 * values are in one stretch, and all the negative values in another - then you'd be adding 
 		 * and subtracting large values. The best case would be alternating positive and negative 
 		 * values but that's harder to do, as you'd have to mess with sumptr additions, above.
 		 */
		for (int i=0; i<=mdist; ++i) {
			sumptr[i]-=rmean[i]*fmean[i]*(clen-i-nonempty[i]); // A bit of numerical instability with this step.
			sumptr[i]+=rmean[i]*sumfdiff[i]+fmean[i]*sumrdiff[i];
		
			// Providing some protection when sd is zero (e.g. empty post-shift vectors, delay is greater than the chromosome length). 
			if (i>ffirst || i>rfirst) {
				sumptr[i]=0;
			} else {
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

}
