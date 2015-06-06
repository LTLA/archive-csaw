#include "csaw.h"
#include <queue>

enum posttype { START, MIDSTART, MIDEND, END };

struct signpost {
	signpost(int p, posttype t, int l, int i) : position(p), type(t), library(l), index(i) {}
	int position;
	posttype type;
	int library, index;
	bool operator> (const signpost& right) const {
		return (position > right.position);
	}
};

SEXP check_bimodality (SEXP all, SEXP regstart, SEXP regend, SEXP priorcount, SEXP invert) try {
	// Setting structures for the data.
    if (!isNewList(all)) { throw std::runtime_error("data on fragments must be contained within a list"); }
    const int nlibs=LENGTH(all);
	std::deque<const int*> left_ptrs(nlibs), right_ptrs(nlibs), strand_ptrs(nlibs);
	std::deque<int> nums(nlibs), indices(nlibs), widths(nlibs);
	std::priority_queue<signpost, std::deque<signpost>, std::greater<signpost> > next;
	
	for (int i=0; i<nlibs; ++i) {
		SEXP current=VECTOR_ELT(all, i);
		if (!isNewList(current) || LENGTH(current)!=4) { 
			throw std::runtime_error("fragment data must be supplied as a list with start, strand and width"); }
		
		for (int j=0; j<3; ++j) {
			SEXP current_col=VECTOR_ELT(current, j);
			if (!isInteger(current_col)) { throw std::runtime_error("fragment data must be in integer format"); }
			const int* ptr=INTEGER(current_col);
			switch (j) {
				case 0: 
					left_ptrs[i]=ptr; 
					nums[i]=LENGTH(current_col);
					break;
				case 1:
					right_ptrs[i]=ptr;
					if (LENGTH(current_col)!=nums[i]) { throw std::runtime_error("length of vectors must be equal"); }
					break;
				case 2:
					strand_ptrs[i]=ptr;
					if (LENGTH(current_col)!=nums[i]) { throw std::runtime_error("length of vectors must be equal"); }
					break;
				default: 
					break;
			}
		}
		
		SEXP curwidth=VECTOR_ELT(current, 3);
		if (!isInteger(curwidth) || LENGTH(curwidth)!=1) { throw std::runtime_error("width must be an integer vector"); }
		widths[i]=asInteger(curwidth);

		// Populating the priority queue.
		if (nums[i]) { next.push(signpost(left_ptrs[i][0]-widths[i]+1, START, i, 0)); }
	}

	// Setting up structures for the regions.
	if (!isInteger(regstart) || !isInteger(regend)) { throw std::runtime_error("region starts or ends should be integer vectors"); }
	const int nregs=LENGTH(regstart);
	if (nregs!=LENGTH(regend)) { throw std::runtime_error("length of region vectors should be equal"); }
	const int * rs_ptr=INTEGER(regstart), 
		      * re_ptr=INTEGER(regend);
	int reg_index=0, next_regstart=-1;
	if (nregs) { 
		next.push(signpost(rs_ptr[0], START, -1, 0));
		next_regstart=rs_ptr[0];
	}

	// Setting up the prior count and inversion flag.
	if (!isReal(priorcount) || LENGTH(priorcount)!=1) { throw std::runtime_error("double-precision scalar required for the prior count"); }
	const double pc=asReal(priorcount);
	if (!isLogical(invert) || LENGTH(invert)!=1) { throw std::runtime_error("inversion flag should be a logical scalar"); }
	const int inv=asLogical(invert);
	
	SEXP output=PROTECT(allocVector(REALSXP, nregs));
try {
	double* optr=REAL(output);
	
	// Odds and ends.
	std::set<int> current_regs;
	int current_position, current_library, current_index;
	posttype current_type;
	int left_forward=0, left_reverse=0, right_forward=0, right_reverse=0;
	double current_score=1; // default when no reads added; prior counts cancel out.

	std::deque<int> new_regs;
	bool modified_stats=false;
	size_t itnr;
	std::set<int>::iterator itcr;

	// Running through the set; stopping when there are no regions, or when everything is processed.
	while (!next.empty()) { 

		// Pulling out all features at this current position.
		current_position=next.top().position;
//		Rprintf("At position %i\n", current_position);
		do {
			current_library=next.top().library;
			current_type=next.top().type;
			current_index=next.top().index;
//			Rprintf("\t\t processing %i: %i -> %i\n", current_type, current_library, current_index);
			next.pop();

			if (current_library<0) {
				// Using negative library values to encode the region index.
				if (current_type==START) {
					current_regs.insert(current_index);
					new_regs.push_back(current_index);
					if ((++reg_index) < nregs) { 
						next.push(signpost(rs_ptr[reg_index], START, -1, reg_index));
						next_regstart=rs_ptr[reg_index];
					} else {
						next_regstart=-1;
					}
					next.push(signpost(re_ptr[current_index]+1, END, -1, current_index));
				} else {
					current_regs.erase(current_index);
				}

			} else {
				// Recording the number of forward/reverse reads to the left or right of the current position.
				modified_stats=true;
				const int& isforward=strand_ptrs[current_library][current_index];
				switch (current_type) { 
					case START:
						if (isforward) { ++right_forward; }
						else { ++right_reverse; }
						next.push(signpost(left_ptrs[current_library][current_index],
									MIDSTART, current_library, current_index));
						break;
					case MIDSTART:
						if (isforward) { ++left_forward; } 
						else { ++left_reverse; }
						next.push(signpost(right_ptrs[current_library][current_index]+1,
									MIDEND, current_library, current_index));
						break;
					case MIDEND:
						if (isforward) { --right_forward; }
						else { --right_reverse; }
						next.push(signpost(right_ptrs[current_library][current_index] + widths[current_library], 
									END, current_library, current_index));
						break;
					case END:
						if (isforward) { --left_forward; }
						else { --left_reverse; }
						break;
					default:
						break;
				}

				// Adding the next thing (skipping intervening reads that don't affect any regions).
				if (current_type==START) {
					int& next_index=indices[current_library];
					while ((++next_index) < nums[current_library]) { 
						if (current_regs.empty() && next_regstart >=0 &&
								next_regstart >= right_ptrs[current_library][next_index] + widths[current_library]) {
 							continue; 
						}
						next.push(signpost(left_ptrs[current_library][next_index] - widths[current_library]+1, 
									START, current_library, next_index));
						break;
					}
				}
			}
		} while (!next.empty() && next.top().position==current_position);

		// Running through and asking if the updated bimodality score is higher or lower.
		if (modified_stats) {
			if (inv) {
				current_score = std::min((double(left_reverse)+pc)/(double(left_forward)+pc), 				
					(double(right_forward)+pc)/(double(right_reverse)+pc));
			} else {
				current_score = std::min((double(left_forward)+pc)/(double(left_reverse)+pc), 				
					(double(right_reverse)+pc)/(double(right_forward)+pc));
			}
		}
//		Rprintf("\t values are Left forward/reverse: %i/%i, right reverse/forward %i/%i\n", left_forward, left_reverse, right_reverse, right_forward);
//		Rprintf("\t score is %.3f\n", current_score);
		if (!new_regs.empty()) { 
			for (itnr=0; itnr<new_regs.size(); ++itnr) { 
				optr[new_regs[itnr]]=current_score; 
			}
			new_regs.clear();
		}
		if (modified_stats) { 
			for (itcr=current_regs.begin(); itcr!=current_regs.end(); ++itcr) {
 			    if (optr[*itcr] < current_score) { optr[*itcr]=current_score; }
			}
			modified_stats=false;
		}
 	   	
		// Quitting if there's no more regions to process.
		if (current_regs.empty() && next_regstart < 0) { break; }
	}

} catch (std::exception& e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception &e) { 
	return mkString(e.what());
}
