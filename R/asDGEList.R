setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })

setMethod("asDGEList", "SummarizedExperiment", function(object, lib.sizes, norm.factors, assay=1, ...) 
# This defines a wrapper function to convert a SummarizedExperiment class
# object into a DGEList object for input into edgeR.
#
# written by Aaron Lun
# created 2 September 2014
# last modified 3 March 2017
{
	all.args <- list(...)
	if (missing(lib.sizes)) { 
		if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals
	}
	all.args$lib.size <- lib.sizes

	if (missing(norm.factors)) {
		if (!is.null(object$norm.factors)) { 
			all.args$norm.factors <- object$norm.factors 
		}
	} else {
		all.args$norm.factors <- norm.factors 
	}

    all.args$counts <- assay(object, assay=assay) 
	y <- do.call(DGEList, all.args)

    offset <- assays(object)$offset 
    if (!is.null(offset)) {
        y <- scaleOffset(y, offset)
    }
	return(y)
})
