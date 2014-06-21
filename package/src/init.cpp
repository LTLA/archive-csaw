#include "csaw.h"

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
	{"collate_exon_data", (DL_FUNC) &collate_exon_data, 4},
	{"annotate_overlaps", (DL_FUNC) &annotate_overlaps, 13},
	{"best_in_cluster", (DL_FUNC) &best_in_cluster, 3},
	{NULL, NULL, 0}
};

void attribute_visible R_init_csaw(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
};

}
