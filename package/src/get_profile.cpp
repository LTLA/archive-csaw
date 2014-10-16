#include "csaw.h"
#include <set>

/* This function scans through the width and pulls out local maxima,
 * along with the surrounding coverage profile.
 */

SEXP get_profile(SEXP track, SEXP aw, SEXP md) try {
	if (!isInteger(track)) { throw std::runtime_error("binned counts should be an integer vector"); }
	const int * tptr=INTEGER(track);
	const int nbins=LENGTH(track);

	if (!isInteger(md) || LENGTH(md)!=1) { throw std::runtime_error("minimum depth should be an integer scalar"); }
	const int depth=asInteger(md);
	if (!isInteger(aw) || LENGTH(aw)!=1) { throw std::runtime_error("width should be an integer scalar"); }
	const int width=asInteger(aw);

	// Filling up the set to identify the local maxima.
	sort_row_index<int> comp(tptr);
	typedef std::multiset<int, sort_row_index<int> > order_set;
	order_set incoming(comp);
	std::deque<order_set::iterator> ordered;
	for (int i=0; i<=width && i<nbins; ++i) { ordered.push_back(incoming.insert(i)); }

	SEXP output=PROTECT(allocVector(VECSXP, 2));
try {
	SET_VECTOR_ELT(output, 0, allocVector(REALSXP, width));
	double* opptr=REAL(VECTOR_ELT(output, 0));
	for (int i=0; i<width; ++i) { opptr[i]=0; }
	SET_VECTOR_ELT(output, 1, allocVector(INTSXP, width));
	int* onptr=INTEGER(VECTOR_ELT(output, 1));
	for (int i=0; i<width; ++i) { onptr[i]=0; }

	/* We start from the left, and we work our way across the chromosome.  We
 	 * only record a profile if the current element is the local maxima (no
 	 * ties allowed, for the record), and if that profile has a read count
 	 * above the specified minimum.
	 */
	int temp, j;
   	order_set::iterator it;
	for (int i=0; i<nbins; ++i) {
		it=incoming.end();
		--it;

		if (*it==i && tptr[i] >= depth) {
			temp=0;
			if (it!=incoming.begin()) {
				--it;
				if (tptr[*it]==tptr[i]) { 
					temp=1;
				}
			}

			if (!temp) {
				temp=0;
				for (j=i+1; j<=i+width && j<nbins; ++j) {
					opptr[temp]+=double(tptr[j])/tptr[i];
					++(onptr[temp]);
					++temp;
				}
				temp=0;
				for (j=i-1; j>=i-width && j>=0; --j) {
					opptr[temp]+=double(tptr[j])/tptr[i];
					++(onptr[temp]);
					++temp;
				}
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
} catch (std::exception& e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

