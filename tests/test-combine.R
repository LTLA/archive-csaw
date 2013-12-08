# This tests the combining power of the combineFDR function.

nativemerge <- function(reg, tol, sign=NULL) {
	n<-length(reg)
	o<-order(reg)
	reg<-reg[o]
	
	increment<-integer(n)
	chrends<-runLength(seqnames(reg))
	increment[c(1, cumsum(chrends[-length(chrends)])+1)]<-1L
	to.next<-c(0L, start(reg[-1])-end(reg[-n])-1)
	increment[to.next > tol]<-1L

	# If a sign is supplied...
	if (!is.null(sign)) { 
		posfc<-sign[o]
		altered.sign<-c(TRUE, posfc[-1]!=posfc[-n])
		increment[altered.sign]<-1L
	}

	merged.ids<-cumsum(increment)
	merged.ids[o]<-merged.ids
	return(merged.ids)
}

comp <- function(reg, tab, tol) {
	ids<-mergeWindows(reg, tol=tol)
	merged.ids<-nativemerge(reg, tol=tol)
	if (!identical(merged.ids, ids$id)) { stop("merging IDs are not identical") }
	out<-combineFDR(ids$id, tab)

	# Aggregating on the merging ID's to get average log-FC, average log-CPM, and Simes.
	almostidentical <- function(x, y, tol=1e-8) { 
		if (length(x)!=length(y)) { return(FALSE) }
		return(all(abs((x-y)/(x+y+1e-6)) < tol))
	}

	ologfc<-aggregate(tab$logFC~merged.ids, FUN=mean, data=NULL)
	if (!almostidentical(ologfc[,2], out$logFC)) { stop("average log-FC's are not identical"); }
	ologcpm<-aggregate(tab$logCPM~merged.ids, FUN=mean, data=NULL)
	if (!almostidentical(ologcpm[,2], out$logCPM)) { stop("average log-CPM's are not identical"); }
	osimes<-aggregate(tab$PValue~merged.ids, FUN=function(x) { min(p.adjust(x, method="BH")) }, data=NULL)
	if (!almostidentical(osimes[,2], out$PValue)) { stop("combined p-values are not identical"); }
	if (!almostidentical(p.adjust(osimes[,2], method="BH"), out$FDR)) { stop("q-values are not identical"); }

	# Checking if we get the same results after reversing the ids (ensures internal re-ordering is active).
	re.o<-length(ids$id):1
	out2<-combineFDR(ids$id[re.o], tab[re.o,])
	if (!almostidentical(out$logFC, out2$logFC) || !almostidentical(out$logCPM, out2$logCPM)
			|| !almostidentical(out$PValue, out2$PValue)) { stop("values not preserved after shuffling"); }

	# Adding some tests if there's multiple log-FC's in 'tab'.
	is.fc<-which(colnames(tab)=="logFC")
	colnames(tab)[is.fc]<-"logFC.1"
	tab$logFC.2<--tab$logFC.1
    out<-combineFDR(ids$id, tab)
	if (!almostidentical(out$logFC.1, -out$logFC.2)) { stop("check failed for multiple log-FC columns") }

	# Also checking with weights for the Simes p-values.
	rand.weights<-runif(length(ids$id), 0, 1)
    out2<-combineFDR(ids$id, tab, weight=rand.weights)
	checker<-split(data.frame(tab$PValue, rand.weights), merged.ids)
	wsimes<-sapply(checker, FUN=function(x) {
		o<-order(x[,1])
		min(x[o,1]*sum(x[,2])/cumsum(x[o,2]))
	})
	if (!almostidentical(wsimes, out2$PValue)	) { stop("weighted combined p-values don't match up") }

	# Checking the reported value of each region.
	ostarts<-aggregate(start(reg)~ merged.ids, FUN=min, data=NULL)
	if (!identical(ostarts[,2], start(ids$region))) { stop("region starts do not match up") }
	oends<-aggregate(end(reg)~merged.ids, FUN=max, data=NULL)
	if (!identical(oends[,2], end(ids$region))) { stop("region ends do not match up") }
	to.use<-!duplicated(merged.ids)
	if (!identical(seqnames(reg)[to.use], seqnames(ids$region))) { stop("region chromosomes do not match up") }
	return(ids$region)
}

generateWindows <- function(chrs, nwin, winsize) {
	allregs<-GRanges()
	for (x in names(chrs)) {
		max.step<-floor(chrs[[x]]/nwin)
		stopifnot(max.step >= 1)
		pos<-cumsum(round(runif(nwin, 1, max.step)))
		suppressWarnings(allregs<-c(allregs, GRanges(x, IRanges(pos, 
			pmin(chrs[[x]], pos+winsize-1L)))))
	}
	total.n<-nwin*length(chrs)
	tab<-data.frame(logFC=runif(total.n, -1, 1), logCPM=runif(total.n, -2, 1),
 		   PValue=rbeta(total.n, 1, 10))
	return(list(region=allregs, table=tab))
}

###################################################################################################

set.seed(2135045)
suppressPackageStartupMessages(require(csaw))

chrs<-c(chrA=1000, chrB=2000, chrC=5000)
blah<-generateWindows(chrs, 200, 1)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 200, 1)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 200, 10)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 100, 50)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 100, 1)
comp(blah$region, blah$table, 10L)

blah<-generateWindows(chrs, 100, 1)
comp(blah$region, blah$table, 20L)

blah<-generateWindows(chrs, 100, 20)
comp(blah$region, blah$table, 10L)

blah<-generateWindows(chrs, 100, 1)
comp(blah$region, blah$table, 20L)

blah<-generateWindows(chrs, 100, 1)
comp(blah$region, blah$table, 20L)

# Trying again with some larger chromosomes.

chrs<-c(chrA=10000, chrB=20000, chrC=50000)
blah<-generateWindows(chrs, 500, 1)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 200, 1)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 500, 10)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 500, 50)
comp(blah$region, blah$table, 1000L)

blah<-generateWindows(chrs, 500, 1)
comp(blah$region, blah$table, 10L)

blah<-generateWindows(chrs, 500, 1)
comp(blah$region, blah$table, 20L)

blah<-generateWindows(chrs, 500, 20)
comp(blah$region, blah$table, 10L)

blah<-generateWindows(chrs, 500, 20)
comp(blah$region, blah$table, 10L)

blah<-generateWindows(chrs, 500, 20)
comp(blah$region, blah$table, 10L)

###################################################################################################
###################################################################################################
# Sticking some tests for merging of fixed-size windows with switched fold changes.

mcomp <- function(tol=1000, ...) {
	stuff<-generateWindows(...)
	posfc<-stuff$tab$logFC>0
	combo<-mergeWindows(stuff$region, tol=tol, sign=posfc)
	merged.ids<-nativemerge(stuff$region, tol=tol, sign=posfc)
	if (!identical(merged.ids, combo$id)) { stop("merging identities failed to match") }

	signs<-split(posfc, merged.ids)
	if (!all(sapply(signs, FUN=function(x) { length(unique(x))==1L }))) {
		stop("not all windows have the same sign") }
	return(combo$region)
}

mcomp(100, chrs=chrs, nwin=200, winsize=1)
mcomp(100, chrs=chrs, nwin=200, winsize=10)
mcomp(100, chrs=chrs, nwin=200, winsize=100)

mcomp(100, chrs=chrs, nwin=500, winsize=1)
mcomp(100, chrs=chrs, nwin=500, winsize=10)
mcomp(100, chrs=chrs, nwin=500, winsize=100)

mcomp(1000, chrs=chrs, nwin=200, winsize=10)
mcomp(1000, chrs=chrs, nwin=500, winsize=10)
mcomp(1000, chrs=chrs, nwin=600, winsize=10)

# Seeing what happens with nested windows of same and opposing sign.

gr <- GRanges("chrA", IRanges(c(1,1,1), c(100, 30, 50))) # should be okay, start point equality
x <- mergeWindows(gr, tol=10, sign=c(TRUE, TRUE, TRUE))
x <- mergeWindows(gr, tol=10, sign=c(TRUE, FALSE, TRUE))

gr <- GRanges("chrA", IRanges(c(100, 20, 40), c(200, 200, 200))) # should be okay, end point equality
x <- mergeWindows(gr, tol=10, sign=c(TRUE, TRUE, TRUE))
x <- mergeWindows(gr, tol=10, sign=c(TRUE, FALSE, TRUE))

gr <- GRanges("chrA", IRanges(c(1, 3, 50), c(200, 100, 80))) 
x <- mergeWindows(gr, tol=10, sign=c(TRUE, TRUE, TRUE)) # should be okay
try(x <- mergeWindows(gr, tol=10, sign=c(TRUE, FALSE, TRUE))) # should fail.
gr2 <- GRanges(c("chrA", "chrB", "chrA"), IRanges(c(1, 3, 50), c(200, 100, 80))) 
x <- mergeWindows(gr2, tol=10, sign=c(TRUE, FALSE, TRUE)) # should be okay again

# In this case, we should see the use of the larger nested region to match to the nest window,
# as they should be exactly 99 bp apart.
gr <- c(gr, GRanges("chrA", IRanges(300, 400)))
x <- mergeWindows(gr, tol=99)
stopifnot(length(unique(x$id))==1L)
x <- mergeWindows(gr, tol=98)
stopifnot(length(unique(x$id))==2L)

###################################################################################################

