getBestTest <- function(ids, tab, mode=c("PValue", "logCPM"), weight=rep(1, length(ids)))
# This uses Holms' method to provide strong control of the FWER within each
# cluster. The idea is that it returns the test with the lowest p-value in the
# cluster. You can then use one test as the representative of the entire cluster,
# which is more specific than Simes (where it's a vague statement of, the DB
# event is somewhere in there). 
#
# written by Aaron Lun
# 17 April 2014
{
	if (!is.integer(ids)) { ids <- as.integer(ids + 0.5) }
	mode <- match.arg(mode) 
	id.order<-order(ids)
	ids<-ids[id.order]
	tab<-tab[id.order,]

	if (mode=="PValue") { 
		if (!is.double(weight)) { weight <- as.double(weight) }
		if (is.null(tab$PValue)) { stop("result table should have one PValue field") }
		weight <- weight[id.order]
			
		# Identifying the minimum P-value.
		out<-.Call("R_best_in_cluster", tab$PValue, ids, weight, PACKAGE="csaw")
		if (is.character(out)) { stop(out) }
		result<-data.frame(best=id.order[out[[2]]], PValue=out[[1]], FDR=p.adjust(out[[1]], method="BH"))

	} else if (mode=="logCPM") {
		if (is.null(tab$logCPM)) { stop("result table should have one logCPM field") }
		weight <- rep(1, length(ids))

		# Identifying the maximum logCPM.
		out<-.Call("R_best_in_cluster", -tab$logCPM, ids, weight, PACKAGE="csaw")
		if (is.character(out)) { stop(out) }
		pval <- tab$PValue[out[[2]]]
		result<-data.frame(best=id.order[out[[2]]], PValue=pval, FDR=p.adjust(pval, method="BH"))
	}
	
	rownames(result)<-ids[c(TRUE, diff(ids)!=0L)]
	return(result)
}

# You can reduce the conservativeness of the Bonferroni method by using windows that 
# don't overlap. You can also just use combineTests with smaller clusters if you're getting
# lost within each cluster.
#    More complex approaches would attempt to estimate the correlation based on genome-wide
# data. I think the best approach would be to compute the 1st order statistic for p-values
# in each cluster of a given size across the dataset. You can then fit the observed statistics
# robustly (i.e., get rid of very low p-values) to a B(1, n) distribution to determine the
# effective number of independent tests 'n'. You can then apply that to the Bonferroni 
# correction for a cluster of that size.

