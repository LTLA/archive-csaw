# This defines some wrapper functions for various statistical methods 
# in edgeR, for the SummarizedExperiment class. The idea is to relieve
# the need for manual specification of the inputs every time; especially
# for the library sizes.
#
# written by Aaron Lun
# created 2 September 2014

setMethod("normalize", "SummarizedExperiment", function(object, lib.sizes=NULL, ...) {
	if (is.null(lib.sizes)) { 
		if (is.null(object$totals)) {
			return(normalizeCounts(counts=assay(object), ...))
		} 
		lib.sizes <- object$totals
	}
	normalizeCounts(counts=assay(object), lib.sizes=lib.sizes, ...)
})

setGeneric("asDGEList", function(object, lib.sizes=NULL, ...) { standardGeneric("asDGEList") })
setMethod("asDGEList", "SummarizedExperiment", function(object, lib.sizes=NULL, ...) {
	if (is.null(lib.sizes)) { lib.sizes <- object$totals }
	DGEList(assay(object), lib.size=lib.sizes, ...)
})
