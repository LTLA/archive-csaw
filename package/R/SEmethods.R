# This defines some wrapper functions for various statistical methods 
# in edgeR, for the SummarizedExperiment class. The idea is to relieve
# the need for manual specification of the inputs every time; especially
# for the library sizes (some care is required to still allow that).
#
# written by Aaron Lun
# created 2 September 2014
# last modified 26 April 2015

setMethod("normalize", "SummarizedExperiment", function(object, ...) {
	opt <- match.call(normalizeCounts, do.call(call, list("normalizeCounts", 
		counts=assay(object), ...)))

	if (is.null(opt$lib.sizes)) { 
		if (is.null(object$totals)) {
			warning("library sizes not found in 'totals'")
		} else {
			opt$lib.sizes <- object$totals
		}
	}

	eval(opt)
})

setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })
setMethod("asDGEList", "SummarizedExperiment", function(object, ...) {
	opt <- match.call(DGEList, do.call(call, list("DGEList", 
		counts=assay(object), ...)))
	
	if (is.null(opt$lib.size)) { 
		if (is.null(object$totals)) { 
			warning("library sizes not found in 'totals'")	
		} else {
			opt$lib.size <- object$totals
		}
	}
	
	eval(opt)
})
