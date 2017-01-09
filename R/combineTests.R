combineTests <- function(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL)
# Computes a combined FDR by assembling their group numbers and computing the
# average log-FC, average log-CPM and Simes' p-value for each cluster. The idea 
# is to test the joint null for each cluster, and then use that to compute the 
# FDR across all clusters. Clusters can be formed by whatever means are deemed 
# necessary (e.g. mergeWindows below, or using peaks).
# 
# written by Aaron Lun
# created 30 July 2013
# last modified 9 January 2017
{
    input <- .check_test_inputs(ids, tab, weight)
    ids <- input$ids
    tab <- input$tab
    groups <- input$groups
    weight <- input$weight
 
	# Running the clustering procedure.
    fc.col <- .parseFCcol(fc.col, tab) 
    is.pval <- .getPValCol(pval.col, tab)
	out <- .Call(cxx_get_cluster_stats, fc.col - 1L, is.pval - 1L, tab, ids, weight, 0.5)
	if (is.character(out)) { stop(out) }

	combined <- data.frame(out[[1]], out[[2]], out[[3]], p.adjust(out[[3]], method="BH"), row.names=groups)
	colnames(combined) <- c("nWindows", 
			sprintf("%s.%s", rep(colnames(tab)[fc.col], each=2), c("up", "down")), 
			colnames(tab)[is.pval], "FDR")

    # Adding direction.
    if (length(fc.col)==1L) {
        labels <- c("mixed", "up", "down")
        combined$direction <- labels[out[[4]] + 1L]
    }
    return(combined)
}

.check_test_inputs <- function(ids, tab, weight) {
    f <- factor(ids)
    all.names <- levels(f)
    ids <- as.integer(f)
    
	if (is.null(weight)) { 
        weight <- rep(1, length(ids)) 
    } else if (!is.double(weight)) { 
        weight <- as.double(weight) 
    }
	stopifnot(length(ids)==nrow(tab))
	stopifnot(length(ids)==length(weight))

    okay.ids <- !is.na(ids)
    if (!all(okay.ids)) { 
        ids <- ids[okay.ids]
        weight <- weight[okay.ids]
        tab <- tab[okay.ids,]
    }
	id.order <- order(ids)
	ids <- ids[id.order]
	tab <- tab[id.order,]
	weight <- weight[id.order]

    if (!all(okay.ids)) { 
        originals <- which(okay.ids)[id.order]
    } else {
        originals <- id.order
    }
    return(list(ids=ids, groups=all.names, tab=tab, weight=weight, original=originals))
}

.parseFCcol <- function(fc.col, tab, multiple=TRUE) {
    if (is.null(fc.col)) { 
        fc.col <- grep("logFC", colnames(tab))
    } else if (is.character(fc.col)) { 
        fc.col <- match(fc.col, colnames(tab)) 
        if (any(is.na(fc.col))) { stop("failed to match logFC column names") }
    }
    if (!multiple && length(fc.col)!=1L) {
        stop("only one column should be specified by 'fc.col'")
    }
    as.integer(fc.col)
}
