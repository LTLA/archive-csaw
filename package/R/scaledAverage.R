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

getWidths <- function(data, len=NULL)
# This computes the effective width of the data in the SummarizedExperiment
# object. This is done by accounting for the effect of read extension; or, for
# paired end data, the median fragment length (from `rescue.ext`, or 
# from `len` if `rescue.ext` is not specified).
#
# written by Aaron Lun
# created 5 November 2014
# last modified 13 February 2015
{
	if (!is.null(len)) { len <- rep(len, length.out=ncol(data)) }
	is.pe <- sapply(data$param, FUN=function(x) { x$pe=="both" })
	frag.len <- integer(ncol(data))

	# Single-end.
	frag.len[!is.pe] <- data$final.ext[!is.pe]
	missing.se <- is.na(frag.len[!is.pe])
	if (any(missing.se)) {
		if (is.null(len)) { 
			stop("need to specify fragment lengths for single-end data")
		}
		frag.len[!is.pe][missing.se] <- len[!is.pe][missing.se]
	}
	
	# Paired-end.
	pe.len <- sapply(data$param, FUN=function(x) { x$rescue.ext })
	not.def <- is.na(pe.len)
	use.pe.len <- is.pe & !not.def
	frag.len[use.pe.len] <- pe.len[use.pe.len]
	use.def.len <- is.pe & not.def
	if (any(use.def.len)) { 
		if (is.null(len)) { 
			stop("need to specify fragment lengths for paired-end data")
		}
		frag.len[use.def.len] <- len[use.def.len]	
	}

	frag.len <- as.integer(mean(frag.len))
	width(rowData(data)) + frag.len - 1L
}

