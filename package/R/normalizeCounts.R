normalizeCounts <- function(counts, lib.sizes=NULL, type=c("scaling", "loess"), weighted=FALSE, ...) 
# This provides a wrapper to perform TMM normalization with non-standard
# library sizes (e.g. due to filtering) and weighting turned off.
# Alternatively, it can do a form a fast loess-like normalization which uses
# the average count as the covariate, rather than the typical A-value-based
# shenanigans. This avoids instability at low counts.
#
# written by Aaron Lun
# created 19 November 2013
# last modified 22 July 2015
{
	if (is.null(lib.sizes)) { 
		warning("library sizes not specified, column sums used instead")
		lib.sizes <- colSums(counts) 
	}

	type <- match.arg(type)
	if (type=="scaling") { 
		y <- DGEList(counts, lib.size=lib.sizes)
		y <- calcNormFactors(y, doWeighting=weighted, ...)
		return(y$samples$norm.factors)

	} else if (type=="loess") { 
		# Scaled corrections squeeze offsets towards relative log-library sizes.
		# Constant value of `cont.cor' would squeeze them towards zero.
		cont.cor <- 0.5
		cont.cor.scaled <- cont.cor * lib.sizes/mean(lib.sizes)
		
		# Using it as a prior.count for abundance ensures linearity with log-counts.
		ab <- aveLogCPM(counts, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))

		offs <- matrix(0, nrow(counts), ncol(counts), byrow=TRUE)
		for (x in seq_len(ncol(counts))) {
			fit <- loessFit(log(counts[,x]+cont.cor.scaled[x]), ab, ...)
			offs[,x] <- fit$fitted 
		}
		offs <- offs-rowMeans(offs)
		return(offs)
	}
}

