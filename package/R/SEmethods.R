# This defines some wrapper functions for various statistical methods 
# in edgeR, for the SummarizedExperiment class. The idea is to relieve
# the need for manual specification of the inputs every time; especially
# for the library sizes.
#
# written by Aaron Lun
# 2 September, 2014

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

setGeneric("param", function(object) { standardGeneric("param") })
setMethod("param", "SummarizedExperiment", function(object) {
	if (any(object$param==0L)) {
		if (!all(object$param==0L)) {
			stop("zero readParam indices mixed in with non-zero indices")
		}
		return(exptData(object)$param)
	} else {
		if (length(unique(object$param))==1L) {
			return(exptData(object)$param[[object$param[1]]])
		} else {
			return(exptData(object)$param[object$param])
		}
	}
})
