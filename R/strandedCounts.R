strandedCounts <- function(bam.files, param=readParam(forward=NULL), regions=NULL, ...) 
# This is a convenience wrapper for strand-specific window counting.  It
# returns a RangedSummarizedExperiment object with up to two sets of counts
# for each window. This is logistically easier than doing it inside
# windowCounts itself.
#
# written by Aaron Lun
# created 9 February 2015
# last modified 2 December 2015
{
	nbam <- length(bam.files)
    if (length(param$forward)) { stop("set forward=NULL in param for strand-specific counting") } 

	if (is.null(regions)) { 
		fdata <- windowCounts(bam.files=bam.files, param=reform(param, forward=TRUE), ...)
		rdata <- windowCounts(bam.files=bam.files, param=reform(param, forward=FALSE), ...)
	} else {
		fdata <- regionCounts(bam.files=bam.files, param=reform(param, forward=TRUE), regions=regions, ...)
		rdata <- regionCounts(bam.files=bam.files, param=reform(param, forward=FALSE), regions=regions, ...)
	}	
	
	# Combining them together.
	combined <- SummarizedExperiment(rbind(assay(fdata), assay(rdata)),
		rowRanges=c(rowRanges(fdata), rowRanges(rdata)),
		colData=colData(fdata), metadata=metadata(fdata))

	o <- order(rowRanges(combined))
	combined <- combined[o]
	combined$totals <- fdata$totals + rdata$totals
	combined$forward.totals <- fdata$totals
	combined$reverse.totals <- rdata$totals
	return(combined)
}
