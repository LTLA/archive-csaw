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
	is.pe <- sapply(paramList(data), FUN=function(x) { x$pe=="both" })
	frag.len <- integer(ncol(data))
	if (is.na(exptData(data)$final.ext)) { 
		frag.len[!is.pe] <- data$ext[!is.pe]
	} else {
		frag.len[!is.pe] <- exptData(data)$final.ext
	}
	
	pe.len <- sapply(paramList(data), FUN=function(x) { x$rescue.ext })
	not.def <- is.na(pe.len)
	use.pe.len <- is.pe & !not.def
	frag.len[use.pe.len] <- pe.len[use.pe.len]
	use.def.len <- is.pe & not.def
	if (any(use.def.len)) { 
		warning("using a median fragment length of 100 bp for PE data")
		frag.len[use.def.len] <- 100L		
	}

	frag.len <- as.integer(mean(frag.len))
	width(rowData(data)) + frag.len - 1L
}

