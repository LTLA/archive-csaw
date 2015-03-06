filterWindows <- function(data, background, type="global", prior.count=2) 
# This is a function for proportion- or background-based filtering of a
# SummarizedExperiment object. For the former, it computes the relative ranks
# that can be used to determine the proportion of highest-abundance windows to
# keep. For the latter, it returns the enrichment term between data and
# background.
#
# written by Aaron Lun
# created 18 February 2015	
# last modified 6 March 2015
{
	type <- match.arg(type, c("global", "local", "control", "proportion"))
	abundances <- scaledAverage(asDGEList(data), scale=1, prior.count=prior.count)

	if (type=="proportion") {
		genome.windows <- .getWindowNum(data)
		relative.rank <- 1 - (rank(abundances) - 1)/genome.windows
		return(list(abundances=abundances, filter=relative.rank))

	} else {
		if (missing(background) && type=="global") {
			filter.stat <- abundances  - .getGlobalBg(data, abundances, prior.count)
			return(list(abundances=abundances, filter=filter.stat))
		} 

		bwidth <- getWidths(background)
		dwidth <- getWidths(data)

		if (type=="global") { 
			.checkLibSizes(data, background)
			relative.width <- median(bwidth)/median(dwidth)
			bg.ab <- scaledAverage(asDGEList(background), scale=relative.width, prior.count=prior.count)
			filter.stat <- abundances - .getGlobalBg(background, bg.ab, prior.count)
			
		} else if (type=="local") {
 		    if (!identical(nrow(data), nrow(background))) { stop("data and background should be of the same length") }	
			.checkLibSizes(data, background)
			relative.width <- (bwidth  - dwidth)/dwidth		
			bg.y <- asDGEList(background)
			bg.y$counts <- bg.y$counts - assay(data)

			# Some protection for negative widths (counts should be zero, so only the prior gets involved in bg.ab).
			subzero <- relative.width <= 0
			if (any(subzero)) { 
				relative.width[subzero] <- 1
				bg.y$counts[subzero,] <- 0L
			}	
			bg.ab <- scaledAverage(bg.y, scale=relative.width, prior.count=prior.count)
			filter.stat <- abundances - bg.ab

		} else {
 		    if (!identical(nrow(data), nrow(background))) { stop("data and background should be of the same length") }	
			relative.width <- bwidth/dwidth
			lib.adjust <- prior.count * mean(background$totals)/mean(data$totals) # Account for library size differences.
			bg.ab <- scaledAverage(asDGEList(background), scale=relative.width, prior.count=lib.adjust)
			filter.stat <- abundances - bg.ab
		}

		return(list(abundances=abundances, back.abundances=bg.ab, filter=filter.stat))
	}
}

.checkLibSizes <- function(data, background) {
	if (!identical(data$totals, background$totals)) { 
		stop("data and background totals should be identical")
	}
	return(NULL)
}

.getWindowNum <- function(data) 
# Get the total number of windows, to account for those not 
# reported in windowCounts (for empty windows/those lost by filter > 1).
{
	spacing <- exptData(data)$spacing
	if (is.null(spacing)) { stop("failed to find spacing for windows") }
	sum(ceiling(seqlengths(rowRanges(data))/spacing)) 
}

.getGlobalBg <- function(data, ab, prior.count)
# Getting the quantile of those windows that were seen, corresponding to 
# the median of all windows in the genome.
{
	prop.seen <- length(ab)/.getWindowNum(data)
 	if (prop.seen > 1) { return(median(ab)) }
	if (prop.seen < 0.5) { return(aveLogCPM(rbind(integer(ncol(data))), lib.size=data$totals, prior.count=prior.count)) }
	quantile(ab, probs=1 - 0.5/prop.seen) 
}
