upweightSummit <- function(ids, summits)
# This takes a set of cluster IDs and a vector of summits (e.g., from
# getBestTest or from findMaxima). It then computes relative weights for each
# cluster ID, with the best entry upweighted by the number of tests.
#
# written by Aaron Lun
# created 26 February 2015
# last modified 14 January 2016
{
    not.okay.ids <- is.na(ids) 
    not.okay.sum <- is.na(summits)
    if (any(not.okay.ids) || any(not.okay.sum)) { 
        final <- numeric(length(ids))
        final[!not.okay.ids] <- Recall(ids[!not.okay.ids], summits[!not.okay.sum])
        return(final)
    }

    if (!is.integer(ids)) { ids <- as.integer(ids) }
	weights <- rep(1, length(ids))
	if (any(ids <= 0L)) { stop("cluster ID vector should have positive integers") }
	freq <- tabulate(ids)
	sum.ids <- ids[summits]
	freq.sum <- tabulate(sum.ids)
	weights[summits] <- freq[sum.ids]/freq.sum[sum.ids]
	return(weights)
}
