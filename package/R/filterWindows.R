filterWindows <- function(data, background, type="global", prior.count=NULL, len=NULL) 
# This is a function for proportion- or background-based filtering of a
# SummarizedExperiment object. For the former, it computes the relative ranks
# that can be used to determine the proportion of highest-abundance windows to
# keep. For the latter, it returns the enrichment term between data and
# background.
#
# written by Aaron Lun
# created 18 February 2015	
# last modified 21 February 2015
{
	type <- match.arg(type, c("global", "local", "control", "proportion"))
	if (is.null(prior.count)) { 
		if (type %in% c("local", "control")) { # To stabilize imprecise background estimates. 
			prior.count <- 5
		} else {
			prior.count <- 2
		}
	}
	abundances <- scaledAverage(asDGEList(data), scale=1, prior.count=prior.count)

	if (type=="proportion") {
		spacing <- exptData(data)$spacing
		if (is.null(spacing)) { stop("proportion-based filtering only works with window/bin counts") }
		genome.windows <- sum(seqlengths(rowData(data))/spacing) # To account for those lost by filter>1 in windowCounts.
		relative.rank <- 1 - (rank(abundances) - 1)/genome.windows
		return(list(abundances=abundances, filter=relative.rank))

	} else {
		if (missing(background) && type=="global") { 
			filter.stat <- abundances  - median(abundances)
			return(list(abundances=abundances, filter=filter.stat))
		} 

		bwidth <- getWidths(background, len=len)
		dwidth <- getWidths(data, len=len)

		if (type=="global") { 
			relative.width <- median(bwidth)/median(dwidth)
			bg.ab <- scaledAverage(asDGEList(background), scale=relative.width, prior.count=prior.count)
			filter.stat <- abundances - median(bg.ab)
			
		} else if (type=="local") {
 		    if (!identical(nrow(data), nrow(background))) { stop("data and background should be of the same length") }	
			relative.width <- (bwidth  - dwidth)/dwidth		
			bg.y <- asDGEList(background)
			bg.y$counts <- bg.y$counts - assay(data)

			# Some protection for negative widths (counts should be zero, so only the prior gets involved in bg.ab).
			subzero <- relative.width <= 0
			if (any(subzero)) { 
				relative.width[subzero] <- 1
				bg.y$counts[sub.zero,] <- 0
			}	
			bg.ab <- scaledAverage(bg.y, scale=relative.width, prior.count=prior.count)
			filter.stat <- abundances - bg.ab

		} else {
 		    if (!identical(nrow(data), nrow(background))) { stop("data and background should be of the same length") }	
			relative.width <- bwidth/dwidth
			bg.ab <- scaledAverage(asDGEList(background), scale=relative.width, prior.count=prior.count)
			filter.stat <- abundances - bg.ab
		}

		return(list(abundances=abundances, back.abundances=bg.ab, filter=filter.stat))
	}
}
