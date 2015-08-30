setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })

setMethod("asDGEList", "RangedSummarizedExperiment", function(object, lib.sizes, ...) 
# This defines a wrapper function to convert a RangedSummarizedExperiment
# object into a DGEList object for input into edgeR.
#
# written by Aaron Lun
# created 2 September 2014
# last modified 29 August 2015
{
	if (missing(lib.sizes)) { 
		if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals
	}
	DGEList(assay(object), lib.size=lib.sizes, ...)
})
