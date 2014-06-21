#ifndef CSAW_H
#define CSAW_H

extern "C" {
#include "R.h"
#include "Rdefines.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

/* From annotator.cpp */

SEXP collate_exon_data (SEXP, SEXP, SEXP, SEXP);

SEXP annotate_overlaps (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
     SEXP, SEXP, SEXP,
	 SEXP, SEXP, SEXP, SEXP); 

/* From best_in_cluster.cpp */

SEXP best_in_cluster(SEXP, SEXP, SEXP);

}

#include <stdexcept>
#include <deque>
#include <algorithm>

template <class T>
struct sort_row_index {
	sort_row_index(const T* p) : ptr(p) {}
	bool operator() (const int& l, const int& r) const { return (ptr[l] < ptr[r]); }
private:
	const T* ptr;
};


#endif

