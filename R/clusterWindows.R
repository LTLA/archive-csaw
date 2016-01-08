clusterWindows <- function(regions, tab, target=0.05, pval.col=NULL, tol, sign=NULL, ..., grid.param=NULL) 
# This does a search for the clusters based on DB windows. 
# It aims to achieve a cluster-level FDR of 'target'.
#
# written by Aaron Lun
# created 8 January 2016
{
    regions <- .toGRanges(regions)
    if (nrow(tab)!=length(regions)) { stop("number of regions is not consistent with entries in 'tab'") }
    pval.col <- .getPValCol(pval.col, tab)
    adjp <- p.adjust(tab[,pval.col], method="BH")

    FUN <- function(sig) { mergeWindows(regions[sig], tol=tol, sign=sign[sig], ...) }
    out <- controlClusterFDR(target=target, adjp=adjp, FUN=function(sig) { FUN(sig)$id }, grid.param=grid.param)
    sig <- adjp <= out$threshold
    clusters <- FUN(sig)

    full.ids <- rep(NA_integer_, nrow(tab))
    full.ids[sig] <- clusters$id
    clusters$id <- full.ids
    clusters$FDR <- out$FDR 
    return(clusters)
}
