scaledAverage <- function(y, scale=1, prior.count=NULL, ...)
# This computes the scaled average abundance, with some finesse to deal with
# the interaction between scaling and the prior count. The `scale` factor
# represents the downscaling factor for the abundances, so the prior count
# has to be upscaled during the calculation itself.
#
# written by Aaron Lun
# 5 November 2014
{
	if (any(scale <= 0)) { stop("scaling factor should be positive") }
	if (is.null(prior.count)) { prior.count <- formals(aveLogCPM.DGEList)$prior.count }
	aveLogCPM(y, prior.count=scale*prior.count, ...) - log2(scale)
}

getWidths <- function(data)
# This computes the effective width of the data in the SummarizedExperiment
# object. This is done by accounting for the effect of read extension; or,
# for paired end data, the median fragment length. For multiple libraries,
# the mean fragment length is used.
#
# written by Aaron Lun
# created 5 November 2014
# last modified 12 December 2014
{
	is.pet <- sapply(exptData(data)$param, FUN=function(x) { x$pet=="both" })
	is.def <- sapply(exptData(data)$param, FUN=function(x) { !is.null(attr(x$ext, "default")) })
	if (any(is.pet & is.def)) { 
		warning("assuming that ext holds median fragment length for PE data")
	}
	frag.len <- mean(sapply(exptData(data)$param, FUN=function(x) { x$ext }))
	width(rowData(data)) + frag.len - 1L
}

