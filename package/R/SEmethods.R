# This defines some wrapper functions for various statistical methods 
# in edgeR, for the RangedSummarizedExperiment class. The idea is to relieve
# the need for manual specification of the inputs every time; especially
# for the library sizes (some care is required to still allow that).
#
# written by Aaron Lun
# created 2 September 2014
# last modified 27 April 2015

setMethod("normalize", "RangedSummarizedExperiment", function(object, lib.sizes, ...) {
	if (missing(lib.sizes)) { 
		if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals 
	}
	normalizeCounts(assay(object), lib.sizes=lib.sizes, ...)
})

setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })
setMethod("asDGEList", "RangedSummarizedExperiment", function(object, lib.sizes, ...) {
	if (missing(lib.sizes)) { 
		if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals
	}
	DGEList(assay(object), lib.size=lib.sizes, ...)
})
