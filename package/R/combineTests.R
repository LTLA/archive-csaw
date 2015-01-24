combineTests <- function(ids, tab, weight=rep(1, length(ids)), pval.col=NULL, other.col=NULL)
# Computes a combined FDR by assembling their group numbers and computing the
# average log-FC, average log-CPM and Simes' p-value for each cluster. The idea 
# is to test the joint null for each cluster, and then use that to compute the 
# FDR across all clusters. Clusters can be formed by whatever means are deemed 
# necessary (e.g. mergeWindows below, or using peaks).
# 
# written by Aaron Lun
# Created 30 July 2013
# Last modified 25 January 2015
{
	if (!is.integer(ids)) { ids <- as.integer(ids+0.5) }
	if (!is.double(weight)) { weight <- as.double(weight) }
	id.order <- order(ids)
	ids <- ids[id.order]
	tab <- tab[id.order,]
	weight  <-  weight[id.order]

	# Saying which columns are to be averaged over.
	if (length(other.col)==0L) { 
		is.fcs <- grep("logFC", colnames(tab))-1L	
		if (!length(is.fcs)) { stop("result table should have at least one logFC field") }
		is.cpm <- which(colnames(tab)=="logCPM")-1L		
		if (length(is.cpm)!=1L) { stop("result table should have one logCPM field") }
		other.col <- c(is.fcs, is.cpm)
	} else {
		other.col <- as.integer(other.col) - 1L
	}

	# Saying which column is the p-value field.
	if (length(pval.col)==0L) { 
		is.pval <- which(colnames(tab)=="PValue")-1L
		if (length(is.pval)!=1L) { stop("result table should have one PValue field") }
	} else if (length(pval.col)!=1L) {
		stop("only one p-value field is possible")
	} else {
		is.pval <- as.integer(pval.col) - 1L
	}
 
	# Running the clustering procedure.
	out <- .Call(cxx_get_cluster_stats, other.col, is.pval, tab, ids, weight)
	if (is.character(out)) { stop(out) }
	combined <- data.frame(out[[1]], out[[2]], p.adjust(out[[2]], method="BH"), row.names=ids[c(TRUE, diff(ids)!=0L)])
	colnames(combined) <- c(colnames(tab)[other.col+1L], colnames(tab)[is.pval+1L], "FDR")
	return(combined)
}

