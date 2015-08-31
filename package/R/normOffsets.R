setGeneric("normOffsets", function(object, ...) { standardGeneric("normOffsets") })

setMethod("normOffsets", "matrix", function(object, lib.sizes=NULL, type=c("scaling", "loess"), weighted=FALSE, ...) 
# This provides a wrapper to perform TMM normalization with non-standard
# library sizes (e.g. due to filtering) and weighting turned off.
# Alternatively, it can do a form a fast loess-like normalization which uses
# the average count as the covariate, rather than the typical A-value-based
# shenanigans. This avoids instability at low object.
#
# written by Aaron Lun
# created 19 November 2013
# last modified 29 August 2015
{
	if (is.null(lib.sizes)) { 
		warning("library sizes not specified, column sums used instead")
		lib.sizes <- colSums(object) 
	}

	type <- match.arg(type)
	if (type=="scaling") { 
		y <- DGEList(object, lib.size=lib.sizes)
		y <- calcNormFactors(y, doWeighting=weighted, ...)
		return(y$samples$norm.factors)

	} else if (type=="loess") { 
		# Scaled corrections squeeze offsets towards relative log-library sizes.
		# Constant value of `cont.cor' would squeeze them towards zero.
		cont.cor <- 0.5
		cont.cor.scaled <- cont.cor * lib.sizes/mean(lib.sizes)
		
		# Using it as a prior.count for abundance ensures linearity with log-object.
		ab <- aveLogCPM(object, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))

		offs <- matrix(0, nrow(object), ncol(object), byrow=TRUE)
		for (x in seq_len(ncol(object))) {
			fit <- loessFit(log(object[,x]+cont.cor.scaled[x]), ab, ...)
			offs[,x] <- fit$fitted 
		}
		offs <- offs-rowMeans(offs)
		return(offs)
	}
})

setMethod("normOffsets", "RangedSummarizedExperiment", function(object, lib.sizes, ...) {
	if (missing(lib.sizes)) { 
		if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals 
	}
	normOffsets(assay(object), lib.sizes=lib.sizes, ...)
})

setMethod("normalize", "RangedSummarizedExperiment", function(object, lib.sizes, ...) {
    .Deprecated(new="normOffsets", old="normalize")
    normOffsets(object, lib.sizes=lib.sizes, ...)
})
