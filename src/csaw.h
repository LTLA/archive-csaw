#ifndef CSAW_H
#define CSAW_H

extern "C" {
#include "R.h"
#include "Rdefines.h"
}
#include <stdexcept>
#include <deque>
#include <algorithm>

struct sort_row_index {
	sort_row_index(const double* p) : ptr(p) {}
	bool operator() (const int& l, const int& r) const { return (ptr[l] < ptr[r]); }
private:
	const double* ptr;
};

#endif
