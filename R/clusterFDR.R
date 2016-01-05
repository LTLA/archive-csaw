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

controlClusterFDR <- function(target, adjp, FUN, ..., maxiter=10)
# Identifies the window-level FDR threshold that is required to 
# control the cluster-level threshold at 'target', given the 
# window-level adjusted p-values and the clustering function FUN.
#
# written by Aaron Lun
# created 5 January 2016
{
    lt <- log(target/(1-target))
    grid.range <- 10

    # Using an iterative grid search, as this tends to be most
    # robust for a discrete and discontinuous function.
    for (it in seq_len(maxiter)) { 
        grid <- seq(lt-grid.range, lt+grid.range, length=21)
        thresholds <- exp(grid)/(exp(grid)+1)
        fdrs <- integer(length(grid))

        for (tx in seq_along(thresholds)) { 
            threshold <- thresholds[tx]
            fdrs[tx] <- clusterFDR(FUN(adjp <= threshold, ...), threshold)
        }

        # Grid contracts at a moderate pace, to provide better resolution.
        chosen <- which.min((fdrs - target)^2)
        lt <- grid[chosen]
        grid.range <- grid.range/2
    }

    return(list(threshold=exp(lt)/(exp(lt)+1), cluster.FDR=fdrs[chosen]))
}

