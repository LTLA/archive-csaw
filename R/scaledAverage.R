scaledAverage <- function(y, scale=1, prior.count=NULL, dispersion=NULL)
# This computes the scaled average abundance, with some finesse to deal with
# the interaction between scaling and the prior count. The `scale` factor
# represents the downscaling factor for the abundances, so the prior count
# has to be upscaled during the calculation itself.
#
# written by Aaron Lun
# created 5 November 2014
# last modified 12 November 2016
{
    if (!is(y, "DGEList")) { stop("'y' should be a DGEList") }
	if (any(scale <= 0)) { stop("scaling factor should be positive") }
	if (is.null(prior.count)) { prior.count <- formals(aveLogCPM.DGEList)$prior.count }
    if (is.null(dispersion)) { 
        dispersion <- y$common.dispersion 
        if (is.null(dispersion)) dispersion <- 0.05
    }

    # Computing the prior to add, scaling it, and then adding it.
    # This way ensures that the offsets are constant regardless of 'scale'.
    empty <- matrix(0, nrow(y$counts), ncol(y$counts))
    ap <- addPriorCount(empty, lib.size=y$samples$lib.size*y$samples$norm.factors, prior.count=prior.count)
    ap$y <- ap$y * scale
    ap$y <- ap$y + y$counts

    # Computing the average abundances in ave-logCPM.
	ave <- mglmOneGroup(y=ap$y, offset=ap$offset, dispersion=dispersion, weights=y$weights) 
    ave <- (ave - log(scale) + log(1e6))/log(2) 
    return(ave)
}

getWidths <- function(data) 
# This computes the effective width of the data in the
# RangedSummarizedExperiment object. This is done by accounting for the effect
# of read extension; or, for paired end data, the median fragment length.
#
# written by Aaron Lun
# created 5 November 2014
# last modified 17 December 2015
{
	flen <- metadata(data)$final.ext

	if (is.na(flen)) { 
		flen <- data$ext
		is.missing <- is.na(flen)
		if (any(is.missing)) { 
			if (is.null(data$rlen) || any(is.na(data$rlen[is.missing]))) { 
				stop("need to specify read lengths in 'data$rlen'")
			}
			flen[is.missing] <- data$rlen[is.missing]
		}
		flen <- as.integer(mean(flen))
	}

	width(rowRanges(data)) + flen - 1L
}

