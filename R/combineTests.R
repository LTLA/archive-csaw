combineTests <- function(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL)
# Computes a combined FDR by assembling their group numbers and computing the
# average log-FC, average log-CPM and Simes' p-value for each cluster. The idea 
# is to test the joint null for each cluster, and then use that to compute the 
# FDR across all clusters. Clusters can be formed by whatever means are deemed 
# necessary (e.g. mergeWindows below, or using peaks).
# 
# written by Aaron Lun
# created 30 July 2013
# last modified 8 January 2016
{
	if (!is.integer(ids)) { ids <- as.integer(ids+0.5) }
	if (is.null(weight)) { weight <- rep(1, length(ids)) }
	else if (!is.double(weight)) { weight <- as.double(weight) }
	stopifnot(length(ids)==nrow(tab))
	stopifnot(length(ids)==length(weight))

	id.order <- order(ids)
	ids <- ids[id.order]
	tab <- tab[id.order,]
	weight <- weight[id.order]

	# Saying which columns have the log-fold change field.
	if (length(fc.col)==0L) { 
		fc.col <- grep("logFC", colnames(tab)) - 1L	
		if (!length(fc.col)) { stop("result table should have at least one logFC field") }
	} else {
		if (is.character(fc.col)) { 
			fc.col <- match(fc.col, colnames(tab)) 
			if (any(is.na(fc.col))) { stop("failed to match logFC column names") }
		}
		fc.col <- as.integer(fc.col) - 1L
	}

	# Saying which column is the p-value field.
    is.pval <- .getPValCol(pval.col, tab) - 1L
 
	# Running the clustering procedure.
	out <- .Call(cxx_get_cluster_stats, fc.col, is.pval, tab, ids, weight, 0.5)
	if (is.character(out)) { stop(out) }
	combined <- data.frame(out[[1]], out[[2]], out[[3]], p.adjust(out[[3]], method="BH"))
 	if (length(ids)) { rownames(combined) <- ids[c(TRUE, diff(ids)!=0L)] } 
	colnames(combined) <- c("nWindows", 
			paste0(rep(colnames(tab)[fc.col+1L], each=2), ".", c("up", "down")), 
			colnames(tab)[is.pval+1L], "FDR")

	return(combined)
}

