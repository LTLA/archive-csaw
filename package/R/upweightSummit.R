upweightSummit <- function(ids, summits)
# This takes a set of cluster IDs and a table produced by getBestTest
# with by.pval=FALSE. It then computes relative weights for each cluster
# ID, with the best entry upweighted by the number of tests.
#
# written by Aaron Lun
# created 26 February 2015
{
	weights <- rep(1, length(ids))
	freq <- tabulate(ids)
	sum.ids <- ids[summits]
	freq.sum <- tabulate(sum.ids)
	weights[summits] <- freq[sum.ids]/freq.sum[sum.ids]
	return(weights)
}
