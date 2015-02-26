consolidateSizes <- function(data.list, result.list, equiweight=TRUE, FUN=NULL, 
    merge.args=list(tol=1000), combine.args=list(), region) 
# Consolidates results for multiple window sizes into a result for the
# genomic region over which those windows are tiled. Returns the combined
# results, as well as ID vectors for cross-referencing and inspection.
#
# written by Aaron Lun
# created 26 February 2015
{
	nset <- length(data.list)
	if (nset!=length(result.list)) { stop("data list must have same length as result list") }
	has.integer <- has.sumexp <- FALSE
	nall <- 0L
	for (x in 1:nset) {
		if (is.integer(data.list[[x]])) { 
			currows <- length(data.list[[x]])
			has.integer <- TRUE
		} else {
			currows <- nrow(data.list[[x]])
			has.sumexp <- TRUE
		}
		ntab <- nrow(result.list[[x]])
		if (currows!=ntab) { stop("corresponding entries of data and result lists must have same number of entries") }
		nall <- nall + ntab
	}
	if (has.integer==has.sumexp) { stop("mixing of integer vectors and SummarizedExperiments is not supported") }

	# Merging windows.	
	if (!has.integer) { 
		all.ranges <- list()
		for (x in 1:nset) { all.ranges[[x]] <- rowData(data.list[[x]]) }
		all.ranges <- do.call(c, all.ranges)
		if (is.null(FUN)) { FUN <- function(x) { do.call(mergeWindows, c(merge.args, regions=x)) } }
		merged <- FUN(all.ranges)
		if (! ("id" %in% names(merged) && "region" %in% names(merged)) ) { stop("FUN must yield a list with elements 'id' and 'region'") }
		if (any(merged$id <= 0L)) { stop("merged IDs must be positive integers") }
	} else {
		if (missing(region)) { 
			warning("region should be supplied when clustering is manually performed") 
			region <- NA
		}
		merged <- list(id=unlist(data.list), region=region)
	}

	# Calculating weights, so each window size (or spacing) has the same contribution to the final outcome.
	if (equiweight) {
		last <- 0L
		rel.weights <- numeric(nall)
		for (x in 1:nset) {
			currows <- nrow(result.list[[x]])
			curdex <- last + 1:currows
			curid <- merged$id[curdex]
			ref.weight <- 1/tabulate(curid)
			rel.weights[curdex] <- ref.weight[curid]
			last <- last + currows
		}
	} else { 
		rel.weights <- rep(1, nall) 
	}

	# Combining statistics.
	tabres <- do.call(rbind, result.list)
	tabcom <- do.call(combineTests, c(list(ids=merged$id, tab=tabres, weight=rel.weights), combine.args))

	# Returning output.
	if (!has.integer) { 
		ids <- list()
		last <- 0L
		for (x in 1:nset) { 
			currows <- nrow(result.list[[x]])
			ids[[x]] <- merged$id[last+1:currows]
			last <- last + currows
		}
		names(ids) <- names(data.list)
	} else { 
		ids <- data.list 
	}

	return(list(id=ids, region=merged$region, table=tabcom))
}
