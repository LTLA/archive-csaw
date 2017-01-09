mixedClusters <- function(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL) 
# Tests for mixed clusters, by performing an IUT on the one-sided combined p-values 
# in each direction for each cluster.
#
# written by Aaron Lun
# created 9 January 2017    
{
    fc.col <- .parseFCcol(fc.col, tab, multiple=FALSE)
    pval.col <- .getPValCol(pval.col, tab)
    pval.colname <- colnames(tab)[pval.col]
    all.p <- .make_one_sided(tab, pval.col=pval.col, fc.col=fc.col)

    # Combining the one-sided p-values.
    up.tab <- tab
    up.tab[,pval.col] <- all.p$up
    up.com <- combineTests(ids, up.tab, weight=weight, pval.col=pval.col, fc.col=fc.col)
    up.com$direction <- NULL
    
    # Repeating in the other direction.
    down.tab <- tab
    down.tab[,pval.col] <- all.p$down
    down.com <- combineTests(ids, down.tab, weight=weight, pval.col=pval.col, fc.col=integer(0))
    
    # Taking the IUT p-value.
    up.com[,pval.colname] <- pmax(up.com[,pval.colname], down.com[,pval.colname])
    up.com$FDR <- p.adjust(up.com[,pval.colname], method="BH")
    return(up.com)
}
