# This defines some wrapper functions for various statistical methods 
# in edgeR, for the SummarizedExperiment class. The idea is to relieve
# the need for manual specification of the inputs every time; especially
# for the library sizes (some care is required to still allow that).
#
# written by Aaron Lun
# created 2 September 2014
# last modified 26 April 2015

setMethod("normalize", "SummarizedExperiment", function(object, ...) {
	opt <- list(...)
	is.libsize <- pmatch(names(opt), "lib.sizes")

	if (all(is.na(is.libsize))) { 
		if (is.null(object$totals)) {
			warning("library sizes not found in 'totals'")
		} else {
			opt$lib.sizes <- object$totals
		}
	}

	opt$counts <- assay(object)
	do.call(normalizeCounts, opt)
})

setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })
setMethod("asDGEList", "SummarizedExperiment", function(object, ...) {
	opt <- list(...)
	is.libsize <- pmatch(names(opt), "lib.size")

	if (all(is.na(is.libsize))) { 
		if (is.null(object$totals)) { 
			warning("library sizes not found in 'totals'") 			
		} else {
			opt$lib.size <- object$totals
		}
	}

	opt$counts <- assay(object)
	do.call(DGEList, opt)
})
