#include "csaw.h"

int split_cluster(const int*, const int*, const int&, const int&, const int&, const int&, int*);

void fillSEXP(SEXP&, const int);

/* We assume that incoming elements are sorted by chr -> start -> end. We then proceed 
 * to aggregate elements by chaining together elements that are less than 'tolerance'
 * apart and, if required, have the same 'sign'. We also split them if the difference
 * between the first start and the last end is greater than 'max_size'.
 */

SEXP merge_windows(SEXP chrs, SEXP start, SEXP end, SEXP sign, SEXP tolerance, SEXP max_size) try {
	if (!isInteger(chrs)) { throw std::runtime_error("chromosomes should be a integer vector"); }
	if (!isInteger(start)) { throw std::runtime_error("start vector should be integer"); }
	if (!isInteger(end)) { throw std::runtime_error("end vector should be integer"); }
	if (!isLogical(sign)) { throw std::runtime_error("sign vector should be logical"); }
	if (!isInteger(tolerance) || LENGTH(tolerance)!=1) { throw std::runtime_error("tolerance should be an integer scalar"); }

	// Setting up pointers.
	const int* cptr=INTEGER(chrs);
	const int* sptr=INTEGER(start);
	const int* eptr=INTEGER(end);
	const int* lptr=LOGICAL(sign);
	const int tol=asInteger(tolerance);
	
	// Checking whether we need to supply a maximum size.
	if (!isInteger(max_size) || LENGTH(max_size) > 1) { throw std::runtime_error("maximum size should be an integer scalar"); }
	const bool limit_size=(LENGTH(max_size)==1);
	const int maxs=(limit_size ? asInteger(max_size) : 0);
	
	// Providing some protection against an input empty list.
	const int n = LENGTH(chrs);
	if (n!=LENGTH(start) || n!=LENGTH(end) || n!=LENGTH(sign)) { throw std::runtime_error("lengths of vectors are not equal"); }
		
	// Proceeding with the merge operation.
	SEXP output=PROTECT(allocVector(VECSXP, 4));
	try {
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, n));
		if (n==0) {
			fillSEXP(output, 0);
			UNPROTECT(1);
			return output;
		}
		int* optr=INTEGER(VECTOR_ELT(output, 0));
		int start_index=0, last_end=*eptr;
		bool diffchr, diffsign;
		int i, ngroups;
		*optr=ngroups=1;

		for (i=1; i<n; ++i) {
			diffchr=(cptr[i]!=cptr[i-1]);
 		   	diffsign=(lptr[i]!=lptr[i-1]);

			if (diffchr 											// Obviously, changing if we're on a different chromosome.
				|| sptr[i]-last_end-1 > tol							// Space between windows, start anew if this is greater than the tolerance.
				|| diffsign 										// Checking if the sign is consistent.
		   	) {
 			    if (limit_size) { ngroups=split_cluster(sptr, eptr, last_end, start_index, i, maxs, optr); } // Splitting the cluster, if desired.
				++ngroups;
				optr[i]=ngroups; 
				start_index=i;
			} else {
				optr[i]=optr[i-1];
			}

			/* Fully nested regions don't have a properly defined interpretation when it comes
 			 * to splitting things by sign. We only support consideration of nested regions where
 			 * either of the boundaries are the same. That can be considered to break the 
 			 * previous stretch if it had opposite sign. Otherwise, the next window would have to
 			 * make the decision to match the parent window or its nested child.
 			 *
 			 * If the start is the same, the window with the earlier end point should be ordered
 			 * first, so that won't pass the first 'if'. This means that it'll only enter with
 			 * a fully nested window. Start and end-point equality might be possible at the ends 
 			 * of chromosomes where trimming enforces sameness, but full nesting should not be observed.
 			 *
 			 * If the nested region has the same sign as the parent, then everything proceeds
 			 * normally i.e. same cluster. We make sure to keep 'last_end' as the parent end in
 			 * such cases. This ensures that the following windows get a change to compute
 			 * distances to the parent end (which should be closer).
 			 */
 		    if (!diffchr && eptr[i] < last_end) {
 			   	if (diffsign) { throw std::runtime_error("fully nested windows of opposite sign are not supported"); } 
			} else { last_end=eptr[i]; }
		}

		// Cleaning up the last cluster, if necessary.
  	  	if (limit_size) { ngroups=split_cluster(sptr, eptr, last_end, start_index, n, maxs, optr); }

		// Now, identifying the chromosome, start and end of each region.
		fillSEXP(output, ngroups);
		int* ocptr=INTEGER(VECTOR_ELT(output, 1));
		int* osptr=INTEGER(VECTOR_ELT(output, 2));
		int* oeptr=INTEGER(VECTOR_ELT(output, 3));
		for (i=0; i<ngroups; ++i) { ocptr[i] = -1; }

		int curgroup;
		for (i=0; i<n; ++i) { 
			curgroup=optr[i]-1;
			if (ocptr[curgroup]<0) { 
				ocptr[curgroup]=cptr[i];
				osptr[curgroup]=sptr[i]; // Sorted by start, remember; only need this once.
				oeptr[curgroup]=eptr[i];
			} else if (oeptr[curgroup] < eptr[i]) { 
				oeptr[curgroup]=eptr[i]; 
			}
		}
	} catch (std::exception& e){
		UNPROTECT(1);
		throw;
	}
	
	UNPROTECT(1);
	return output;		
} catch (std::exception& e) {
	return mkString(e.what());
}

int split_cluster(const int* starts, const int* ends, const int& actual_end, const int& xs, const int& xe, const int& width, int* output) { 
	double full_width=actual_end-starts[xs]+1;
	if (full_width <= width) { return output[xs]; }
	int mult=int(std::ceil(full_width/width));

	/* There can only be `mult` subclusters. At the worst, `cur_diff`
	   will be equal to `actual_end - starts[xs]`. Division by `subwidth` will
	   give an expression of `(actual_end - starts[xs]) / [ (actual_end - starts[xs] + 1) / mult ]`.
	   This will always be less than `mult`, so flooring it will give `mult-1`,
	   i.e., the last index of `instantiated`.
	 */
	std::deque<int> instantiated(mult, 0);
	int output_index=output[xs];
	int i=0;

	// Allocating windows into subclusters, based on their midpoints.
	double subwidth=full_width/mult, cur_diff;
	for (i=xs; i<xe; ++i) {
		cur_diff = double(starts[i]+ends[i])*0.5 - starts[xs];
		output[i] = int(cur_diff/subwidth);	
		if (!instantiated[output[i]]) { instantiated[output[i]] = 1; }
	}

	/* Allocating output indices to the subclusters. This approach avoids
	   situations where you get nonconsecutive cluster indices, e.g., when
	   `tol` is greater than the maximum width.	 
	 */
	for (i=0; i<mult; ++i) { 
		if (!instantiated[i]) { continue; }
		instantiated[i]=output_index;
		++output_index;
	}

	// Assigning indices back to the output vector.
	for (i=xs; i<xe; ++i) { output[i]=instantiated[output[i]]; }

	// Returning the last group index that was used.
	return output_index-1;	
}

void fillSEXP(SEXP& output, const int ngroups) {
	SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ngroups));
	SET_VECTOR_ELT(output, 2, allocVector(INTSXP, ngroups));
	SET_VECTOR_ELT(output, 3, allocVector(INTSXP, ngroups));
	return;
}
