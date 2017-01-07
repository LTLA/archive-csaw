empiricalFDR <- function(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL, neg.down=TRUE) 
# Converts two-tailed p-values to one-tailed p-values, combines them 
# and takes the number of rejections in the "wrong" direction as an 
# estimate of the number of false positives.
#
# written by Aaron Lun
# created 7 December 2017
{
    if (length(fc.col)==0L) { 
        fc.col <- grep("logFC", colnames(tab))
    }
    if (length(fc.col)!=1L) {
        stop("only one column should be specified by 'fc.col'")
    }
    cur.fc <- tab[,fc.col]
    if (neg.down) {
        wrong.dir <- cur.fc < 0
    } else {
        wrong.dir <- cur.fc > 0
    }

    # Converting to one-sided p-values.
    pval.col <- .getPValCol(pval.col, tab)
    new.p <- tab[,pval.col]/2
    new.p[wrong.dir] <- 1 - new.p[wrong.dir]

    # Combining p-values in each direction.
    right.tab <- tab
    right.tab[,pval.col] <- new.p
    right.com <- combineTests(ids, right.tab, weight=weight, pval.col=pval.col, fc.col=fc.col)

    wrong.tab <- tab
    wrong.tab[,pval.col] <- 1-new.p
    wrong.com <- combineTests(ids, wrong.tab, weight=weight, pval.col=pval.col, fc.col=fc.col)

    # Computing empirical FDR (enforcing monotonicity).
    right.comp <- right.com$PValue
    empirical <- findInterval(right.comp, sort(wrong.com$PValue))/seq_along(right.comp)
    o <- order(right.comp, decreasing=TRUE)
    empirical[o] <- cummin(empirical[o])

    out <- data.frame(PValue=right.com$PValue, FDR=pmin(1, empirical))
    rownames(out) <- rownames(right.com)
    return(out)
}


