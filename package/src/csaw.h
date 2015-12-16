#ifndef CSAW_H
#define CSAW_H

#include <set>
#include <stdexcept>
#include <deque>
#include <algorithm>
#include <cmath>

template <class T>
struct sort_row_index {
	sort_row_index(const T* p) : ptr(p) {}
	bool operator() (const int& l, const int& r) const { 
		if (ptr[l]==ptr[r]) { return (l < r); }
		else { return (ptr[l] < ptr[r]); }
	}
private:
	const T* ptr;
};

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

extern "C" {

/* annotator.cpp */

SEXP collate_exon_data (SEXP, SEXP, SEXP, SEXP);

SEXP annotate_overlaps (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
     SEXP, SEXP, SEXP,
	 SEXP, SEXP, SEXP, SEXP); 

/* best_in_cluster.cpp */

SEXP best_in_cluster(SEXP, SEXP, SEXP);

/* get_cluster_stats.cpp */

SEXP get_cluster_stats(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* merge_windows.cpp */

SEXP merge_windows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* correlate_reads.cpp */

SEXP correlate_reads(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* get_rle_counts.cpp */

SEXP get_rle_counts(SEXP, SEXP, SEXP, SEXP, SEXP);

/* get_profile.cpp */

SEXP get_profile(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* find_maxima.cpp */

SEXP find_maxima(SEXP, SEXP, SEXP, SEXP, SEXP);

/* check_bimodality.cpp */

SEXP check_bimodality(SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif

