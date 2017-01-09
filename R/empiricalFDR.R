empiricalFDR <- function(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL, neg.down=TRUE) 
# Converts two-tailed p-values to one-tailed p-values, combines them 
# and takes the number of rejections in the "wrong" direction as an 
# estimate of the number of false positives.
#
# written by Aaron Lun
# created 7 January 2017
# last modified 9 January 2017
{
    fc.col <- .parseFCcol(fc.col, tab, multiple=FALSE)
    pval.col <- .getPValCol(pval.col, tab)
    pval.colname <- colnames(tab)[pval.col]
    
    all.p <- .make_one_sided(tab, pval.col=pval.col, fc.col=fc.col)
    if (neg.down) { 
        all.p <- list(right=all.p$up, wrong=all.p$down)
    } else {
        all.p <- list(right=all.p$down, wrong=all.p$up)
    }

    # Combining one-sided p-values.
    right.tab <- tab
    right.tab[,pval.col] <- all.p$right
    right.com <- combineTests(ids, right.tab, weight=weight, pval.col=pval.col, fc.col=fc.col)
    right.com$direction <- NULL
    
    # Repeating in the other direction.
    wrong.tab <- tab
    wrong.tab[,pval.col] <- all.p$wrong
    wrong.com <- combineTests(ids, wrong.tab, weight=weight, pval.col=pval.col, fc.col=integer(0))

    # Computing empirical FDR.
    right.comp <- right.com[,pval.colname]
    o <- order(right.comp)
    right.comp <- right.comp[o]
    empirical <- findInterval(right.comp, sort(wrong.com[,pval.colname]))/seq_along(right.comp)
    
    # Enforcing monotonicity and other characteristics.
    empirical <- pmin(1, empirical)
    empirical <- rev(cummin(rev(empirical)))
    empirical[o] <- empirical
    right.com$FDR <- empirical
    return(right.com)
}

.make_one_sided <- function(tab, pval.col, fc.col) {
    cur.fc <- tab[,fc.col]
    going.up <- cur.fc > 0
    pval <- tab[,pval.col]
    
    # Calculating each set fresh, to avoid numeric imprcesion from repeated "1-" operations
    up.p <- pval/2
    up.p[!going.up] <- 1 - up.p[!going.up]
    down.p <- pval/2
    down.p[going.up] <- 1 - down.p[going.up]
    return(list(up=up.p, down=down.p))
}



