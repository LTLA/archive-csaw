#include "csaw.h"
#include <queue>

typedef std::multiset<int, sort_row_index<double> > order_set;
struct compare_iterators {
	compare_iterators(const int* ptr) : endptr(ptr) {}
	bool operator() (const order_set::iterator& left, const order_set::iterator& right) const {
		if (endptr[*left]==endptr[*right]) { return (*left > *right); } // Opposite, as we want the smallest value to be treated as the largest in the queue.
		return (endptr[*left] > endptr[*right]);
	}
private:
	const int* endptr;
};

/* This function scans through the track and pulls out local maxima. */

SEXP find_maxima(SEXP chrs, SEXP starts, SEXP ends, SEXP metric, SEXP range) try {
	if (!isInteger(chrs) || !isInteger(starts) || !isInteger(ends)) { throw std::runtime_error("chr, start and end vectors must be integer"); }
	const int nlen=LENGTH(starts);
	if (!isReal(metric)) { throw std::runtime_error("metric must be a double-precision vector"); }
	if (LENGTH(chrs)!=nlen || LENGTH(metric) != nlen || LENGTH(ends)!=nlen) { throw std::runtime_error("vectors must be of equal length"); }

	// Pulling out scalars and pointers.	
	if (!isInteger(range) || LENGTH(range)!=1) { throw std::runtime_error("range should be an integer scalar"); }
	const int maxrange=asInteger(range);
	if (maxrange <= 0) { throw std::runtime_error("range should be a positive integer"); }
	const int* sptr=INTEGER(starts), *eptr=INTEGER(ends), *cptr=INTEGER(chrs);
	const double* mptr=REAL(metric);

	const sort_row_index<double> comp(mptr);
	order_set incoming(comp);
	const compare_iterators comp2(eptr);
	typedef std::priority_queue<order_set::iterator, std::deque<order_set::iterator>, compare_iterators> pqueue;
	pqueue first_to_leave(comp2);
	first_to_leave.push(incoming.insert(0));

	SEXP output=PROTECT(allocVector(LGLSXP, nlen));
try {
	if (nlen==0) { 
		UNPROTECT(1);
		return output;
	}
	int* optr=LOGICAL(output);

	// Assuming we're sorted by sptr.
	int right_edge=1, is_max, right_copy;
	double cur_max=mptr[0], max_right;
	order_set::iterator it;
	for (int i=0; i<nlen; ++i) {
		if (i) {
			if (cptr[i] < cptr[i-1]) { throw std::runtime_error("regions must be sorted by chromosome"); }
			else if (cptr[i]==cptr[i-1]) { 
				if (sptr[i] < sptr[i-1]) { throw std::runtime_error("regions on the same chromosome must be sorted by start position"); }
				if (eptr[i] < eptr[i-1]) { throw std::runtime_error("nested regions are not supported"); }
			}
		}

		// Getting rid of regions running off the end.
		while (!first_to_leave.empty() && (cptr[i] != cptr[*(first_to_leave.top())] || sptr[i] - eptr[*(first_to_leave.top())] > maxrange)) {
			incoming.erase(first_to_leave.top());
			first_to_leave.pop();
		}

		// Finding the maximum (records first tie, if there are ties).
		max_right=R_NaReal;
		while (right_edge < nlen && cptr[i]==cptr[right_edge] && sptr[right_edge] - eptr[i] <= maxrange) {
			if (ISNA(max_right) || max_right < mptr[right_edge]) {
				is_max=right_edge;
				max_right=mptr[right_edge];				
			}
			++right_edge;
		}

		if (!ISNA(max_right)) { 	
			/* Checking if the new maximum is greater than the current maximum, in which 
			 * case we clear out everything because the new maximum supercedes anything 
			 * already in the set.
			 */	
			if (cur_max < max_right) {
				incoming.clear();
				first_to_leave = pqueue(comp2);
			}
		
			/* Only adding those after the maximum, and which are also reverse
	 		 * cumulative maxima themselves (ties allowed). There's no point having a 
			 * window which is sandwiched by two maxima, as it'll never be its own maxima.
			 */
			right_copy=right_edge-1;
			max_right=mptr[right_copy];
			while (right_copy >= is_max) { 
				if (mptr[right_copy] >= max_right) { 
					first_to_leave.push(incoming.insert(right_copy));
					max_right=mptr[right_copy];
				}
				--right_copy;
			}
		}

//		Rprintf("Current: %i (%i) %.3f\n", i, right_edge, mptr[i]);
//		for (order_set::const_iterator itx=incoming.begin(); itx!=incoming.end(); ++itx) { 
//			Rprintf("\t%i %.3f\n", *itx, mptr[*itx]);
//		}

		// Checking if we're currently the max (some allowance for exactly tied values).
		optr[i]=0;
		if (incoming.size()) { 
			it=incoming.end();
			--it;
			cur_max=mptr[*it];
			do {
				if (*it==i) { 
					optr[i]=1;
					break;
				}
				if (it==incoming.begin()) { break; } 
				--it;
			} while (mptr[*it]==cur_max);
		} else {
			throw std::runtime_error("empty set during maxima detection");
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

