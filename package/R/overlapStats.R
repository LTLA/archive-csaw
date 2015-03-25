.overlapStats <- function(olap, tab, relation.weights=NULL, test.weights=NULL, type=c("combine", "best"), ...) {
	region.dex <- queryHits(olap)
	win.dex <- subjectHits(olap)

	# Setting up weights.
	if (is.null(relation.weights)) { 
		if (!is.null(test.weights)) {
			relation.weights <- test.weights[win.dex]
		}
	}

	type <- match.arg(type)
	if (type=="combine") { 
		output <- combineTests(region.dex, tab[win.dex,], weight=relation.weights, ...)
	} else { 
		output <- getBestTest(region.dex, tab[win.dex,], weight=relation.weights, ...)
		output$best <- win.dex[output$best]
	}

	# Filling in empties with NA's.
	nregions <- queryLength(olap)
	expand.vec <- rep(NA, nregions)
	row.dex <- as.integer(rownames(output))
	if (any(row.dex <= 0L | row.dex > nregions)) { stop("cluster IDs are not within [1, nregions]") }
	expand.vec[row.dex] <- 1:nrow(output)
	output <- output[expand.vec,]
	rownames(output) <- NULL

	return(output)
}

combineOverlaps <- function(olap, tab, relation.weights=NULL, test.weights=NULL, ...) 
# Wrapper around combineTests for Hits from findOverlaps,
# when windows are overlapped with regions.
#
# written by Aaron Lun
# created 25 March 2015
{
	.overlapStats(olap, tab, relation.weights=relation.weights, 
		test.weights=test.weights, type="combine", ...)
}

getBestOverlaps <- function(olap, tab, relation.weights=NULL, test.weights=NULL, ...) 
# Wrapper around getBestTest for Hits from findOverlaps,
# when windows are overlapped with reigons.
#
# written by Aaron Lun
# created 25 March 2015
{
	.overlapStats(olap, tab, relation.weights=relation.weights, 
		test.weights=test.weights, type="best", ...)
}

summitOverlaps <- function(olap, region.best, test.summits) 
# Wrapper around upweightSummits for Hits from findOverlaps.
#
# written by Aaron Lun
# created 25 March 2015
{
	region.dex <- queryHits(olap)
	win.dex <- subjectHits(olap)

	if (!missing(region.best)) { 
		summit.dex <- region.best[region.dex]
		summits <- !is.na(summit.dex) & win.dex==summit.dex
	} else if (!missing(test.summits)) { 
		if (is.integer(test.summits)) { 
			out <- logical(max(win.dex, test.summits))
			out[test.summits] <- TRUE
			test.summits <- out
		}
		summits <- test.summits[win.dex]
	} else {
		stop("either region.best or test.summits must be specified")
	}

	upweightSummit(region.dex, summits)
}
