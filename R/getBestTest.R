getBestTest <- function(ids, tab, by.pval=TRUE, weight=NULL, pval.col=NULL, cpm.col=NULL)
# This uses Holms' method to provide strong control of the FWER within each
# cluster. The idea is that it returns the test with the lowest p-value in the
# cluster. You can then use one test as the representative of the entire cluster,
# which is more specific than Simes (where it's a vague statement of, the DB
# event is somewhere in there). 
#
# written by Aaron Lun
# created 17 April 2014
# last modified 8 January 2017
{
    input <- .check_test_inputs(ids, tab, weight)
    ids <- input$ids
    tab <- input$tab
    groups <- input$groups
    weight <- input$weight

    pval.col <- .getPValCol(pval.col, tab)
	if (by.pval) { 
		# Identifying the minimum P-value, and Bonferroni-correcting it.
		out <- .Call(cxx_best_in_cluster, tab[,pval.col], ids, weight)
		if (is.character(out)) { stop(out) }
		pval <- out[[1]]
		best <- out[[2]]

	} else {
		if (is.null(cpm.col)) { cpm.col <- "logCPM" }
		if (length(cpm.col)!=1L) { 
			stop("absent or multiple logCPM columns are not supported")
		} 

		# Identifying the window with the maximum logCPM.
		out <- .Call(cxx_best_in_cluster, -tab[,cpm.col], ids, weight)
		if (is.character(out)) { stop(out) }
		best <- out[[2]]
		pval <- tab[best, pval.col]
	}
	
	subtab <- tab[best,]
	subtab[,pval.col] <- pval
	result <- data.frame(best=input$original[best], subtab, FDR=p.adjust(pval, method="BH"), row.names=groups)
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

.getPValCol <- function(pval.col, tab) {
    if (length(pval.col)>1L) { 
        stop("multiple p-value columns are not supported")
    }
	if (is.null(pval.col)) { 
		pval.col <- "PValue"
    }
    if (is.character(pval.col)) {
        pval.col <- which(colnames(tab)==pval.col)
    } else {
        pval.col <- as.integer(pval.col) # coerce to integer, just in case.
    }
    if (length(pval.col)==0) { 
        stop("failed to find any p-value column")
    }
    return(pval.col)
}
