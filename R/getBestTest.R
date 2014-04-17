getBestTest <- function(ids, tab, weight=rep(1, length(ids)))
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
	if (!is.double(weight)) { weight <- as.double(weight) }
	id.order<-order(ids)
	ids<-ids[id.order]
	tab<-tab[id.order,]
	if (is.null(tab$PValue)) { stop("result table should have one PValue field") }
	weight <- weight[id.order]

	# Running the clustering procedure.
	out<-.Call("R_best_in_cluster", tab$PValue, ids, weight, PACKAGE="csaw")
	if (is.character(out)) { stop(out) }
	combined<-data.frame(best=id.order[out[[2]]+1L], PValue=out[[1]], FDR=p.adjust(out[[1]], method="BH"))
	rownames(combined)<-ids[c(TRUE, diff(ids)!=0L)]
	return(combined)
}
