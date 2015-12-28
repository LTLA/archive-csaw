upweightSummit <- function(ids, summits)
# This takes a set of cluster IDs and a vector of summits (e.g., from
# getBestTest or from findMaxima). It then computes relative weights for each
# cluster ID, with the best entry upweighted by the number of tests.
#
# written by Aaron Lun
# created 26 February 2015
# last modified 25 March 2015
{
	if (!is.integer(ids)) { ids <- as.integer(ids + 0.5) }
	weights <- rep(1, length(ids))
	if (any(ids <= 0L)) { stop("cluster ID vector should have positive integers") }
	freq <- tabulate(ids)
	sum.ids <- ids[summits]
	freq.sum <- tabulate(sum.ids)
	weights[summits] <- freq[sum.ids]/freq.sum[sum.ids]
	return(weights)
}
