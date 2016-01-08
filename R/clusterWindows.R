clusterWindows <- function(regions, tab, target=0.05, pval.col=NULL, tol, sign=NULL, ..., weight=NULL, grid.param=NULL) 
# This does a search for the clusters based on DB windows. 
# It aims to achieve a cluster-level FDR of 'target'.
#
# written by Aaron Lun
# created 8 January 2016
{
    regions <- .toGRanges(regions)
    if (nrow(tab)!=length(regions)) { stop("number of regions is not consistent with entries in 'tab'") }
    if (missing(tol)) {
        tol <- 100
        warning("'tol' for 'mergeWindows' set to a default of 100 bp")
    }

    # Computing a frequency-weighted adjusted p-value.
    pval.col <- .getPValCol(pval.col, tab)
    if (is.null(weight)) { weight <- rep(1, nrow(tab)) }
    adjp <- .weightedFDR(tab[,pval.col], weight)

    # Controlling the cluster-level FDR
    FUN <- function(sig) { mergeWindows(regions[sig], tol=tol, sign=sign[sig], ...) }
    out <- controlClusterFDR(target=target, adjp=adjp, FUN=function(sig) { FUN(sig)$id }, 
                             weight=weight, grid.param=grid.param)
    sig <- adjp <= out$threshold
    clusters <- FUN(sig)

    # Cleaning up
    full.ids <- rep(NA_integer_, nrow(tab))
    full.ids[sig] <- clusters$id
    clusters$id <- full.ids
    clusters$FDR <- out$FDR 
    return(clusters)
}

.weightedFDR <- function(p, w) {
    if (length(p)!=length(w)) { stop("weight and p-value vector are not of same length") }
    o <- order(p)
    p <- p[o]
    w <- w[o]
    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*p/cumsum(w))))
    pmin(adjp, 1)
}

consolidateClusters <- function(data.list, result.list, equiweight=TRUE, ...) 
# This does the same as clusterWindows, but for results from many different analyses
# (ostensibly with different window sizes).
#
# written by Aaron Lun
# created 8 January 2016
{
    nset <- length(data.list)
    set.it.vec <- seq_len(nset)
    if (nset!=length(result.list)) { stop("data list must have same length as result list") }
    
    for (x in set.it.vec) {
        data.list[[x]] <- .toGRanges(data.list[[x]])
        currows <- length(data.list[[x]])
        ntab <- nrow(result.list[[x]])
        if (currows!=ntab) { stop("corresponding entries of data and result lists must have same number of entries") }
    }
    
    # Merging everyone together.
    all.data <- do.call(c, data.list)
    all.result <- do.call(rbind, result.list)
    groupings <- rep(seq_along(data.list), lengths(data.list)) 
    
    # Computing weights based on number of windows; each analysis contributes same effective number of tests.
    if (equiweight) { 
        weights <- rep(1/lengths(data.list), lengths(data.list))
    } else {
        weights <- NULL
    }   

    out <- clusterWindows(all.data, all.result, weight=weights, ...)
    out$id <- split(out$id, groupings)
    names(out$id) <- names(data.list)
    return(out)
}

