combineTests <- function(ids, tab, weight=rep(1, length(ids)))
# Computes a combined FDR by assembling their group numbers and computing the
# average log-FC, average log-CPM and Simes' p-value for each cluster. The idea 
# is to test the joint null for each cluster, and then use that to compute the 
# FDR across all clusters. Clusters can be formed by whatever means are deemed 
# necessary (e.g. mergeWindows below, or using peaks).
# 
# written by Aaron Lun
# 30 July 2013
{
	if (!is.integer(ids)) { ids<-as.integer(ids+0.5) }
	if (!is.double(weight)) { weight<-as.double(weight) }
	id.order<-order(ids)
	ids<-ids[id.order]
	tab<-tab[id.order,]

	# Saying which columns are what.
	is.fcs<-grep("logFC", colnames(tab))-1L	
	if (!length(is.fcs)) { stop("result table should have at least one logFC field") }
	is.cpm<-which(colnames(tab)=="logCPM")-1L		
	if (length(is.cpm)!=1L) { stop("result table should have one logCPM field") }
	is.pval<-which(colnames(tab)=="PValue")-1L
	if (length(is.pval)!=1L) { stop("result table should have one PValue field") }
 
	# Running the clustering procedure.
	out<-.Call("R_get_cluster_stats", is.fcs, is.cpm, is.pval, tab, ids, weight, PACKAGE="csaw")
	if (is.character(out)) { stop(out) }
	combined<-data.frame(out[[1]], logCPM=out[[2]], PValue=out[[3]], FDR=p.adjust(out[[3]], method="BH"), rownames=ids[c(TRUE, diff(ids)!=0L)])
	colnames(combined)[1:length(is.fcs)] <- colnames(tab)[is.fcs+1L]
	return(combined)
}

