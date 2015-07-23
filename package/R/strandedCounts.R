strandedCounts <- function(bam.files, param=readParam(forward=NULL), regions=NULL, ...) 
# This is a convenience wrapper for strand-specific window counting.  It
# returns a RangedSummarizedExperiment object with up to two sets of counts
# for each window. This is logistically easier than doing it inside
# windowCounts itself.
#
# written by Aaron Lun
# created 9 February 2015
# last modified 22 July 2015
{
	nbam <- length(bam.files)
	plist <- .makeParamList(nbam, param)
	for (i in seq_len(nbam)) {
		if (length(plist[[i]]$forward)) { stop("set forward=NULL in param for strand-specific counting") } 
	}

	if (is.null(regions)) { 
		fdata <- windowCounts(bam.files=bam.files, param=reformList(plist, forward=TRUE), ...)
		rdata <- windowCounts(bam.files=bam.files, param=reformList(plist, forward=FALSE), ...)
	} else {
		fdata <- regionCounts(bam.files=bam.files, param=reformList(plist, forward=TRUE), 
			regions=regions, ...)
		rdata <- regionCounts(bam.files=bam.files, param=reformList(plist, forward=FALSE), 
			regions=regions, ...)
	}	
	
	# Combining them together.
	combined <- SummarizedExperiment(rbind(assay(fdata), assay(rdata)),
		rowRanges=c(rowRanges(fdata), rowRanges(rdata)),
		colData=colData(fdata), metadata=metadata(fdata))

	o <- GenomicRanges::order(rowRanges(combined))
	combined <- combined[o]
	combined$totals <- fdata$totals + rdata$totals
	combined$forward.totals <- fdata$totals
	combined$reverse.totals <- rdata$totals
	return(combined)
}
