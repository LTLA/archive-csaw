#include "csaw.h"

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
	{"collate_exon_data", (DL_FUNC) &collate_exon_data, 4},
	{"annotate_overlaps", (DL_FUNC) &annotate_overlaps, 13},
	{"best_in_cluster", (DL_FUNC) &best_in_cluster, 3},
	{"get_cluster_stats", (DL_FUNC) &get_cluster_stats, 6},
	{"merge_windows", (DL_FUNC) &merge_windows, 6},
	{"correlate_reads", (DL_FUNC) &correlate_reads, 6},
	{"get_rle_counts", (DL_FUNC) &get_rle_counts, 5}, 
	{"get_profile", (DL_FUNC) &get_profile, 6}, 
	{"find_maxima", (DL_FUNC) &find_maxima, 5}, 
	{"check_bimodality", (DL_FUNC) &check_bimodality, 5}, 
	{NULL, NULL, 0}
};

void attribute_visible R_init_csaw(DllInfo *dll) {
	R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

}
