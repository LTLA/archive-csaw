scaledAverage <- function(y, scale=1, prior.count=NULL, ...)
# This computes the scaled average abundance, with some finesse to deal with
# the interaction between scaling and the prior count. The `scale` factor
# represents the downscaling factor for the abundances, so the prior count
# has to be upscaled during the calculation itself.
#
# written by Aaron Lun
# created 5 November 2014
{
	if (any(scale <= 0)) { stop("scaling factor should be positive") }
	if (is.null(prior.count)) { prior.count <- formals(aveLogCPM.DGEList)$prior.count }
	aveLogCPM(y, prior.count=scale*prior.count, ...) - log2(scale)
}

getWidths <- function(data) 
# This computes the effective width of the data in the
# RangedSummarizedExperiment object. This is done by accounting for the effect
# of read extension; or, for paired end data, the median fragment length.
#
# written by Aaron Lun
# created 5 November 2014
# last modified 14 May 2015
{
	flen <- metadata(data)$final.ext

	if (is.na(flen)) { 
		flen <- data$ext
		is.missing <- is.na(flen)
		if (any(is.missing)) { 
			if (is.null(data$rlen)) { 
				stop("need to specify read lengths in 'data$rlen'")
			}
			flen[is.missing] <- data$rlen[is.missing]
		}
		flen <- as.integer(mean(flen))
	}

	width(rowRanges(data)) + flen - 1L
}

