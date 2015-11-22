setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })

setMethod("asDGEList", "SummarizedExperiment0", function(object, lib.sizes, norm.factors, ...) 
# This defines a wrapper function to convert a SummarizedExperiment class
# object into a DGEList object for input into edgeR.
#
# written by Aaron Lun
# created 2 September 2014
# last modified 22 November 2015
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

    if ("counts" %in% assayNames(object)) {
        all.args$counts <- assays(object)$counts
    } else {
        all.args$counts <- assay(object) # Just using the first, if unnamed.
    }

	y <- do.call(DGEList, all.args)
	y$offset <- assays(object)$offset # just NULL, if it's not there.
	return(y)
})
