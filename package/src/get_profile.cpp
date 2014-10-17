#include "csaw.h"
#include <set>

#ifdef DEBUG
#include <map>
#endif

/* This function scans through the track and pulls out local maxima. */

SEXP get_profile(SEXP starts, SEXP ends, SEXP total_pts, SEXP res, SEXP aw, SEXP md) try {
	
	// Height of stinginess; just calling it from within C++. Avoid need for type checking.
	SEXP output, track=PROTECT(get_rle_counts(starts, ends, total_pts, res, ScalarLogical(1)));
	try{
		if (isString(track)) { 
			throw std::runtime_error(CHAR(STRING_ELT(track, 0)));
		}
		const int * tptr=INTEGER(track);
		const int nbins=LENGTH(track);

		if (!isInteger(md) || LENGTH(md)!=1) { throw std::runtime_error("minimum depth should be an integer scalar"); }
		const int depth=asInteger(md);
		if (!isInteger(aw) || LENGTH(aw)!=1) { throw std::runtime_error("width should be an integer scalar"); }
		const int width=asInteger(aw);
		const int binsize=asInteger(res);

		// Filling up the set to identify the local maxima.
		sort_row_index<int> comp(tptr);
		typedef std::multiset<int, sort_row_index<int> > order_set;
		order_set incoming(comp);
		std::deque<order_set::iterator> ordered;
		for (int i=0; i<=width && i<nbins; ++i) { ordered.push_back(incoming.insert(i)); }

		/* We start from the left, and we work our way across the chromosome. We
	 	 * determine if the current element is the local maxima and if that maxima
	 	 * has a read count above the specified minimum.
	 	 */
		int temp, j;
   		order_set::iterator it;
		std::deque<int> collected, collected_depth;

		for (int i=0; i<nbins; ++i) {
			it=incoming.end();
			--it;

			if (*it==i && tptr[i] >= depth) {
				temp=0;
				if (it!=incoming.begin()) {
					--it;

					// No ties allowed.
					if (tptr[*it]==tptr[i]) { temp=1; }
				}
				if (!temp) { 
					collected.push_back(i); 
					collected_depth.push_back(tptr[i]);
				}
			}

			// Updating the structures for the next position.
			if (i >= width) { 
				incoming.erase(ordered.front());
				ordered.pop_front();
			}
			temp=i+width+1;
			if (temp < nbins) {
				ordered.push_back(incoming.insert(temp));
			}
		}

		/************ Collating the coverage profile ************/

		const int * sptr = INTEGER(starts);
		const int * eptr= INTEGER(ends);
		const int nreads=LENGTH(starts);
		const int totals=asInteger(total_pts);

#ifdef DEBUG
		std::map<int, std::deque<int> > debugger;
		for (int i=0; i<collected.size(); ++i) { debugger[collected[i]].resize(width+1); }
#endif

		output=PROTECT(allocVector(VECSXP, 2));
		try {
			SET_VECTOR_ELT(output, 0, allocVector(REALSXP, width));
			double * ocptr=REAL(VECTOR_ELT(output, 0)) - 1;
			for (int i=1; i<=width; ++i) { ocptr[i]=0; }

			int left, right, diff, counter, refcounter;
			std::deque<int>::const_iterator colstart=collected.begin(), colend=collected.end();
			const int colsize=collected.size();

			for (int i=0; i<nreads; ++i) {
				left=(sptr[i] < 2 ? 0 : int((sptr[i]-2)/binsize)+1);
				right=int((eptr[i]-1)/binsize);

				/* Find relevant maxima, and calculate the location in the
 				 * profile in which these reads lie. We ignore maxima that the
 				 * read overlaps (i.e., negative or zero diff values), because
 				 * we're trying to determine when and how many new reads are
 				 * included when the window is extended.
				 */
				counter = refcounter = std::lower_bound(colstart, colend, left)-colstart;
				
				if (counter!=0) {
					do { 
						--counter;
						diff = left - collected[counter];
						if (diff!=0) { 
							if (diff > width) { break; }
							ocptr[diff] += 1/double(collected_depth[counter]);

#ifdef DEBUG
							int& tango = debugger[collected[counter]][diff];
							if (tango%2==0) { ++tango; }
#endif
						}
					} while (counter > 0);
				}
				
				counter=refcounter;
				while (counter!=colsize) { 
					diff = collected[counter] -  right;
					if (diff > 0) { 
						if (diff > width) { break; }
						ocptr[diff] += 1/double(collected_depth[counter]);

#ifdef DEBUG
						int& tango = debugger[collected[counter]][diff];
						if (tango < 2) { tango+=2; }
#endif
					}
					++counter;
				}
			}

#ifdef DEBUG
			// Printing the debugger:
			for (std::map<int, std::deque<int> >::iterator itx=debugger.begin(); itx!=debugger.end(); ++itx) {
				Rprintf("%i: ", itx->first);
				for (int m=0; m<=width; ++m) { 
					Rprintf("%i, ", (itx->second)[m]); 
				}
				Rprintf("\n");
			}
#endif

			// Getting the contribution of each maxima to the coverage profile (to account for end-of-chromosome truncations).
			SET_VECTOR_ELT(output, 1, allocVector(INTSXP, width));
			int * onptr=INTEGER(VECTOR_ELT(output, 1)) - 1;
			for (int i=1; i<=width; ++i) { onptr[i]=0; }

			for (int i=0; i<colsize; ++i) {
				if (collected[i] < width) {
					if (collected[i]!=0) { ++(onptr[collected[i]]); }
				} else {
					++(onptr[width]);
				}
				
				diff=totals - collected[i] - 1;
				if (diff < width) {
					if (diff!=0) { ++(onptr[diff]); }
				} else {
					++(onptr[width]);
				}
			}	

			for (int i=width-1; i>0; --i) { onptr[i]+=onptr[i+1]; }

		} catch (std::exception& e) {
			UNPROTECT(1);
			throw;
		}
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}

	UNPROTECT(2);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

