clusterFDR <- function(ids, threshold, weight=NULL)
# This computes an informal estimate of the cluster-level FDR,
# given the cluster IDs for all significant windows. The idea
# is to allow clustering of significant windows to explicitly
# identify the differential subinterval in complex regions.
#
# written by Aaron Lun
# created 13 April 2015
# last modified 8 January 2016
{
    ids <- as.integer(ids)
    o <- order(ids)
	ids <- ids[o]

    if (is.null(weight)) { 
        weight <- rep(1, length(ids)) 
    } else {
        if (length(weight)!=length(ids)) { 
            stop("lengths of 'weight' and 'ids' must be the same")
        }
        weight <- as.double(weight)
        weight <- weight[o]
    }

	num.fp <- sum(weight) * threshold
	cluster.sizes <- .Call(cxx_get_cluster_weight, ids, weight) 
	num.fp.cluster <- sum(cumsum(sort(cluster.sizes)) <= num.fp)

    if (length(cluster.sizes)) { 
    	return(num.fp.cluster/length(cluster.sizes))
    } else {
        return(0)
    }
}

controlClusterFDR <- function(target, adjp, FUN, ..., weight=NULL, grid.param=NULL)
# Identifies the window-level FDR threshold that is required to 
# control the cluster-level threshold at 'target', given the 
# window-level adjusted p-values and the clustering function FUN.
#
# written by Aaron Lun
# created 5 January 2016
# last modified 8 January 2016
{
    if (is.list(grid.param)) { 
        ref.args <- c("length", "iter", "range", "scale")
        m <- pmatch(names(grid.param), ref.args)
        if (any(is.na(m))) { stop("invalid arguments specified in 'grid.param'") }
        names(grid.param) <- ref.args[m]
    }

    grid.range <- grid.param$range
    if (is.null(grid.range)) { grid.range <- 20 }
    grid.range <- grid.range/2
    grid.length <- grid.param$length
    if (is.null(grid.length)) { grid.length <- 21 }
    iter <- grid.param$iter
    if (is.null(iter)) { iter <- 5 }
    scale <- grid.param$scale
    if (is.null(scale)) { scale <- 4 }
    
    lt <- log(target/(1-target))
    if (is.null(weight)) { weight <- rep(1, length(adjp)) } 

    # Using an iterative grid search, as this tends to be most
    # robust for a discrete and discontinuous function.
    for (it in seq_len(iter)) { 
        grid <- seq(lt-grid.range, lt+grid.range, length=grid.length)
        thresholds <- exp(grid)/(exp(grid)+1)

        # Doesn't make sense to have window-level FDR > cluster-level FDR. 
        # We'll be conservative and only search grid points that are lower.
        keep <- thresholds <= target
        grid <- grid[keep]
        thresholds <- thresholds[keep]

        fdrs <- integer(length(grid))
        for (tx in seq_along(thresholds)) { 
            threshold <- thresholds[tx]
            is.sig <- adjp <= threshold
            fdrs[tx] <- clusterFDR(FUN(is.sig, ...), threshold, weight=weight[is.sig])
        }

        # Picking the largest threshold that yields an FDR closest to the target.
        # Avoids searching values that are too small when there are many ties.        
        chosen <- length(grid) - which.min(rev(abs(fdrs-target))) + 1L
        lt <- grid[chosen]

        # Grid contracts at a moderate pace, to provide better resolution.
        grid.range <- grid.range/scale
    }

    return(list(threshold=exp(lt)/(exp(lt)+1), FDR=fdrs[chosen]))
}

