clusterFDR <- function(ids, threshold)
# This computes an informal estimate of the cluster-level FDR,
# given the cluster IDs for all significant windows. The idea
# is to allow clustering of significant windows to explicitly
# identify the differential subinterval in complex regions.
#
# written by Aaron Lun
# created 13 April 2015
{
	ids <- sort(ids)
	num.fp <- length(ids) * threshold
	cluster.sizes <- rle(ids)$lengths
	num.fp.cluster <- sum(cumsum(sort(cluster.sizes)) <= num.fp) 
	return(num.fp.cluster/length(cluster.sizes))
}
