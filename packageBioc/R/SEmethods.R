# This defines some wrapper functions for various statistical methods 
# in edgeR, for the SummarizedExperiment class. The idea is to relieve
# the need for manual specification of the inputs every time; especially
# for the library sizes.
#
# written by Aaron Lun
# 2 September, 2014

setGeneric("average", function(object, ...) { standardGeneric("average") })
setMethod("average", "SummarizedExperiment", function(object, ...) {
	if (is.null(object$totals)) { 
		warning("no library sizes found in totals for column data")
	}
	aveLogCPM(assay(object), lib.size=object$totals, ...)
})

setMethod("normalize", "SummarizedExperiment", function(object, ...) {
	if (is.null(object$totals)) { 
		return(normalizeCounts(counts=assay(object), ...))
	}
	normalizeCounts(counts=assay(object), lib.sizes=object$totals, ...)
})

setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })
setMethod("asDGEList", "SummarizedExperiment", function(object, ...) {
	DGEList(assay(object), lib.size=object$totals, ...)
})
