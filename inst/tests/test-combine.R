# This tests the combining power of the combineTests function.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))

comp <- function(total.n, n.clusters, weights=NULL) {
	merged.ids <- sample(n.clusters, total.n, replace=TRUE)
	tab <- data.frame(logFC=runif(total.n, -1, 1), logCPM=runif(total.n, -2, 1),
		PValue=rbeta(total.n, 1, 10))
	if (is.null(weights)) { weights <- rep(1, length.out=total.n) }
	out <- combineTests(merged.ids, tab, weight=weights)

	# Checking numbers.
	if (!identical(as.integer(table(merged.ids)), out$nWindows)) { stop("total numbers of windows are not identical") }
	reference <- integer(nrow(out))
	test <- table(merged.ids[tab$logFC > 0.5])
	reference[match(names(test), rownames(out))] <- as.integer(test)
	if (!identical(reference, out$logFC.up)) { stop("number of up windows is not identical") }
	reference <- integer(nrow(out))
	test <- table(merged.ids[tab$logFC < -0.5])
	reference[match(names(test), rownames(out))] <- as.integer(test)
	if (!identical(reference, out$logFC.down)) { stop("number of up windows is not identical") }
   	
	# Checking Simes.
	checker <- split(data.frame(PValue=tab$PValue, weight=weights), merged.ids)
	osimes<-sapply(checker, FUN=function(x) {
		o<-order(x$PValue)
		min(x$PValue[o]/cumsum(x$weight[o])) * sum(x$weight)
	})
	almostidentical <- function(x, y, tol=1e-8) { 
		if (length(x)!=length(y)) { return(FALSE) }
		return(all(abs((x-y)/(x+y+1e-6)) < tol))
	}
	if (!almostidentical(osimes, out$PValue)) { stop("combined p-values are not identical") }
	if (!almostidentical(p.adjust(osimes, method="BH"), out$FDR)) { stop("q-values are not identical") }

	# Checking the rownames.
	if (!identical(rownames(out), as.character(sort(unique(merged.ids))))) { stop("row names are not matched") }

	# Checking if we get the same results after reversing the ids (ensures internal re-ordering is active).
	re.o <- total.n:1
	if (is.null(weights)) { 
		out2<-combineTests(merged.ids[re.o], tab[re.o,])
	} else {
		out2<-combineTests(merged.ids[re.o], tab[re.o,], weight=weights[re.o])
	}
	if (!almostidentical(out$logFC, out2$logFC) || !almostidentical(out$logCPM, out2$logCPM)
		|| !almostidentical(out$PValue, out2$PValue)) { stop("values not preserved after shuffling") }

	# Adding some tests if there's multiple log-FC's in 'tab'.
	is.fc<-which(colnames(tab)=="logFC")
	colnames(tab)[is.fc]<-"logFC.1"
	tab$logFC.2<--tab$logFC.1
    out<-combineTests(merged.ids, tab, weight=weights)
	if (!identical(out$logFC.1.up, out$logFC.2.down)) { stop("check failed for multiple log-FC columns") }
	if (!identical(out$logFC.1.down, out$logFC.2.up)) { stop("check failed for multiple log-FC columns") }

	return(head(out))
}

###################################################################################################

set.seed(2135045)
suppressWarnings(suppressPackageStartupMessages(require(csaw)))
options(width=100)

comp(20, 5)
comp(20, 10)
comp(20, 20)
comp(20, 5, weights=runif(20))
comp(20, 10, weights=runif(20))
comp(20, 20, weights=runif(20))

comp(100, 50)
comp(100, 100)
comp(100, 200)
comp(100, 50, weights=runif(100))
comp(100, 100, weights=runif(100))
comp(100, 200, weights=runif(100))

comp(1000, 50)
comp(1000, 100)
comp(1000, 200)
comp(1000, 50, weights=runif(1000))
comp(1000, 100, weights=runif(1000))
comp(1000, 200, weights=runif(1000))

###################################################################################################
# Checking for sane behaviour when no IDs are supplied.

combineTests(integer(0), data.frame(PValue=numeric(0), logCPM=numeric(0), logFC=numeric(0)), weight=numeric(0))

###################################################################################################
# End.

