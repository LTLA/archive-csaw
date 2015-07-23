consolidateSizes <- function(data.list, result.list, equiweight=TRUE, 
    merge.args=list(tol=1000), combine.args=list(), region=NULL, overlap.args=list()) 
# Consolidates results for multiple window sizes into a result for the
# genomic region over which those windows are tiled. Returns the combined
# results, as well as ID vectors for cross-referencing and inspection.
#
# written by Aaron Lun
# created 26 February 2015
# last modified 22 July 2015
{
	nset <- length(data.list)
	set.it.vec <- seq_len(nset)
	if (nset!=length(result.list)) { stop("data list must have same length as result list") }
	for (x in set.it.vec) {
		currows <- nrow(data.list[[x]])
		ntab <- nrow(result.list[[x]])
		if (currows!=ntab) { stop("corresponding entries of data and result lists must have same number of entries") }
	}

	# Merging windows, or finding overlaps.
	if (is.null(region)) { 
		all.ranges <- list()
		for (x in set.it.vec) { all.ranges[[x]] <- rowRanges(data.list[[x]]) }
		all.ranges <- do.call(c, all.ranges)
		merged <- do.call(mergeWindows, c(merge.args, regions=all.ranges)) 

		# Formatting for nice output.
		final.ids <- list()
		last <- 0L
		for (x in set.it.vec) { 
			currows <- nrow(result.list[[x]])
			final.ids[[x]] <- merged$id[last+seq_len(currows)]
			last <- last + currows
		}
		names(final.ids) <- names(data.list)
	} else {
		all.ranges <- list()
		final.ids <- list()
		for (x in set.it.vec) {
			olap <- do.call(findOverlaps, c(query=region, subject=rowRanges(data.list[[x]]), overlap.args))
			final.ids[[x]] <- olap
			all.ranges[[x]] <- queryHits(olap)
			result.list[[x]] <- result.list[[x]][subjectHits(olap),]
		}
		merged <- list(id=unlist(all.ranges), region=region)
		names(final.ids) <- names(data.list)
	}

	# Calculating weights, so each window size (or spacing) has the same contribution to the final outcome.
	if (equiweight) {
		last <- 0L
		rel.weights <- list()
		for (x in set.it.vec) {
			currows <- nrow(result.list[[x]])
			curid <- merged$id[last + seq_len(currows)]
			rel.weights[[x]] <- (1/tabulate(curid))[curid]	
			last <- last + currows		
		}
		rel.weights <- unlist(rel.weights)
	} else { 
		rel.weights <- NULL
	}

	# Combining statistics.
	tabres <- do.call(rbind, result.list)
	tabcom <- do.call(combineTests, c(list(ids=merged$id, tab=tabres, weight=rel.weights), combine.args))
	if (!is.null(region)) { tabcom <- .expandNA(tabcom, length(region)) }

	return(list(id=final.ids, region=merged$region, table=tabcom))
}
