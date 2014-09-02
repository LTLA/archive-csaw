normalizeCounts <- function(counts, lib.sizes, type=c("scaling", "loess"), weighted=FALSE, dispersion=0.05, ...) 
# This provides a wrapper to perform TMM normalization with non-standard
# library sizes (e.g. due to filtering) and weighting turned off.
# Alternatively, it can do a form a fast loess-like normalization which uses
# the average count as the covariate, rather than the typical A-value-based
# shenanigans. This avoids instability at low counts.
#
# written by Aaron Lun
# 19 November, 2013
# modified 1 September, 2014
{
	if (missing(lib.sizes)) { 
		lib.sizes <- colSums(counts)
		warning("library sizes not specified, column sums used instead")
	}

	type <- match.arg(type)
	if (type=="scaling") { 
		y<-DGEList(counts, lib.size=lib.sizes)
		y<-calcNormFactors(y, doWeighting=weighted, ...)
		return(y$samples$norm.factors)
	} else if (type=="loess") { 
	    ab <- aveLogCPM(counts, lib.size=lib.sizes, dispersion=dispersion)
		offs <- matrix(0, nrow(counts), ncol(counts), byrow=TRUE)
		for (x in 1:ncol(counts)) { offs[,x]<-loessFit(log(counts[,x]+0.5), ab, ...)$fitted }
		offs<-offs-rowMeans(offs)
		return(offs)
	}
}

# Writing a wrapper for SummarizedExperiment inputs.
setMethod("normalize", "SummarizedExperiment", function(object, ...) {
	if (is.null(object$totals)) { 
		return(normalizeCounts(counts=assay(object), ...))
	}
	normalizeCounts(counts=assay(object), lib.sizes=object$totals, ...)
})

