getBestTest <- function(ids, tab, by.pval=TRUE, weight=NULL, pval.col=NULL, cpm.col=NULL)
# This uses Holms' method to provide strong control of the FWER within each
# cluster. The idea is that it returns the test with the lowest p-value in the
# cluster. You can then use one test as the representative of the entire cluster,
# which is more specific than Simes (where it's a vague statement of, the DB
# event is somewhere in there). 
#
# written by Aaron Lun
# created 17 April 2014
# last modified 25 March 2015
{
	if (!is.integer(ids)) { ids <- as.integer(ids + 0.5) }
	stopifnot(length(ids)==nrow(tab))
	id.order <- order(ids)
	ids <- ids[id.order]
	tab <- tab[id.order,]

	# Checking what's what.
	if (length(pval.col)==0L) { 
		pval.col <- which(colnames(tab)=="PValue")
		if (length(pval.col)==0L) { stop("result table should have one PValue field") }
	} else if (length(pval.col)>1L) { 
		stop("multiple p-value columns are not supported")
	} else { 
		if (is.character(pval.col)) {
			pval.col <- match(pval.col, colnames(tab)) 
			if (any(is.na(pval.col))) { stop("failed to match p-value column names") }
		}
		pval.col <- as.integer(pval.col) 
	}

	if (by.pval) { 
		# Identifying the minimum P-value, and Bonferroni-correcting it.
		if (is.null(weight)) { weight <- rep(1, length(ids)) } 
		else if (!is.double(weight)) { weight <- as.double(weight) }
		stopifnot(length(ids)==length(weight))

		weight <- weight[id.order]
		out <- .Call(cxx_best_in_cluster, tab[,pval.col], ids, weight)
		if (is.character(out)) { stop(out) }
		pval <- out[[1]]
		best <- out[[2]]

	} else {
		if (length(cpm.col)==0L) {
			cpm.col <- which(colnames(tab)=="logCPM")
			if (length(cpm.col)==0L) { stop("result table should have one logCPM field") }
		} else if (length(cpm.col)>1L) { 
			stop("multiple logCPM columns are not supported")
		} else {
			if (is.character(cpm.col)) {
				cpm.col <- match(cpm.col, colnames(tab)) 
				if (any(is.na(cpm.col))) { stop("failed to match CPM column names") }
			}
			cpm.col <- as.integer(cpm.col)
		}
		weight <- rep(1, length(ids))

		# Identifying the window with the maximum logCPM.
		out <- .Call(cxx_best_in_cluster, -tab[,cpm.col], ids, weight)
		if (is.character(out)) { stop(out) }
		best <- out[[2]]
		pval <- tab[best, pval.col]
	}
	
	subtab <- tab[best,]
	subtab[,pval.col] <- pval
	result <- data.frame(best=id.order[best], subtab, FDR=stats::p.adjust(pval, method="BH"))
	if (length(ids)) { rownames(result) <- ids[c(TRUE, diff(ids)!=0L)] }

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

