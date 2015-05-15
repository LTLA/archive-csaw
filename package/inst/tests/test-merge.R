# This tests the mergeWindows function, separately from the combineTests function.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))

nativemerge <- function(reg, tol, sign=NULL) {
	n<-length(reg)
	o<-order(reg)
	reg<-reg[o]
	
	increment<-integer(n)
	last.end <- end(reg)
	by.chr <- split(1:n, as.character(seqnames(reg)))
	for (x in by.chr) {
		increment[x[1]] <- 1L
		last.end[x] <- cummax(last.end[x])
	}
	to.next <- c(0L, start(reg)[-1]-last.end[-n]-1)
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

basecomp <- function(tol=100, ...) {
	reg <- generateWindows(...)
	merged.ids <- nativemerge(reg, tol)
	ids <- mergeWindows(reg, tol=tol)
	if (!identical(merged.ids, ids$id)) { stop("clustering IDs do not match up") }

	# Checking the reported value of each region.
	ostarts<-aggregate(start(reg)~ merged.ids, FUN=min, data=NULL)
	if (!identical(ostarts[,2], start(ids$region))) { stop("region starts do not match up") }
	oends<-aggregate(end(reg)~merged.ids, FUN=max, data=NULL)
	if (!identical(oends[,2], end(ids$region))) { stop("region ends do not match up") }
	if (!identical(seqnames(reg), seqnames(ids$region[merged.ids]))) { stop("region chromosomes do not match up") }

	return(ids$region)
}

source("simsam.R")

###############################################################################################

set.seed(123213)

chrs <- c(chrA=10000, chrB=5000, chrC=2000)
basecomp(chrs=chrs, nwin=50, winsize=1)
basecomp(chrs=chrs, nwin=50, winsize=10)
basecomp(chrs=chrs, nwin=100, winsize=10)

basecomp(chrs=chrs, nwin=50, winsize=1, tol=50)
basecomp(chrs=chrs, nwin=50, winsize=10, tol=50)
basecomp(chrs=chrs, nwin=100, winsize=10, tol=50)

basecomp(chrs=chrs, nwin=100, winsize=runif(100, 1, 50), tol=10)
basecomp(chrs=chrs, nwin=200, winsize=runif(200, 5, 50), tol=10)
basecomp(chrs=chrs, nwin=500, winsize=runif(500, 5, 50), tol=10)

basecomp(chrs=chrs, nwin=100, winsize=runif(100, 5, 50), tol=5)
basecomp(chrs=chrs, nwin=200, winsize=runif(200, 5, 50), tol=5)
basecomp(chrs=chrs, nwin=500, winsize=runif(500, 5, 50), tol=5)

###################################################################################################
# Sticking some tests for merging of fixed-size windows with switched fold changes.

mcomp <- function(tol=100, ...) {
	stuff<-generateWindows(...)
	posfc <- rbinom(length(stuff), 1, 0.5)==1L
	combo<-mergeWindows(stuff, tol=tol, sign=posfc)
	merged.ids<-nativemerge(stuff, tol=tol, sign=posfc)
	if (!identical(merged.ids, combo$id)) { stop("merging identities failed to match for variable sign") }

	signs<-split(posfc, merged.ids)
	if (!all(sapply(signs, FUN=function(x) { length(unique(x))==1L }))) {
		stop("not all windows have the same sign") }

	# Checking that there's actually something being tested, here.
	ref <- nativemerge(stuff, tol=tol)
	if (identical(ref, merged.ids)) { stop("simulations don't have any splitting by sign"); }
	return(combo$region)
}

mcomp(100, chrs=chrs, nwin=200, winsize=1)
mcomp(100, chrs=chrs, nwin=200, winsize=10)
mcomp(100, chrs=chrs, nwin=200, winsize=100)

mcomp(100, chrs=chrs, nwin=500, winsize=1)
mcomp(100, chrs=chrs, nwin=500, winsize=10)
mcomp(100, chrs=chrs, nwin=500, winsize=100)

mcomp(100, chrs=chrs, nwin=200, winsize=10)
mcomp(100, chrs=chrs, nwin=500, winsize=10)
mcomp(100, chrs=chrs, nwin=600, winsize=10)

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
# Testing the maximum limit.

maxcomp <- function(tol=100, maxd=200, ...) {
	stuff<-generateWindows(...)
	combo<-mergeWindows(stuff, tol=tol, max.width=maxd)
	merged.ids<-nativemerge(stuff, tol=tol)

	# Parsing the merged.ids, and splitting them.
	gunk <- split(1:length(merged.ids), merged.ids)
	allstarts <- start(stuff)
	allends <- end(stuff)
	last.id <- 1L
	has.split <- FALSE

	for (x in names(gunk)) {
		ix <- as.integer(x)
		all.dexes <- gunk[[x]]
		curstarts <- allstarts[all.dexes]
		curends <- allends[all.dexes]
		
		full.width <- max(curends)-min(curstarts)+1L
		mult <- ceiling(full.width/maxd)
		if (mult >= 2L) { has.split <- TRUE }
		subwidth <- full.width/mult

		mid.dist <- (curstarts + curends)*0.5 - min(curstarts)
		subcluster <- floor(mid.dist / subwidth)
		new.ids <- as.integer(factor(subcluster)) 

		merged.ids[all.dexes] <- last.id + new.ids - 1L
		last.id <- last.id + max(new.ids)
	}

	if (!has.split) { stop("simulation doesn't involve any split regions") }
	if (!identical(merged.ids, combo$id)) { stop("merging identities failed to match with maximum limit") }
	return(combo$region)
}

maxcomp(50, maxd=500, chrs=chrs, nwin=200, winsize=1)
maxcomp(50, maxd=500, chrs=chrs, nwin=200, winsize=10)
maxcomp(10, maxd=500, chrs=chrs, nwin=200, winsize=100)

maxcomp(50, chrs=chrs, nwin=500, winsize=1)
maxcomp(50, chrs=chrs, nwin=500, winsize=10)
maxcomp(10, chrs=chrs, nwin=500, winsize=100)

maxcomp(10, maxd=500, chrs=chrs, nwin=200, winsize=runif(200, 1, 100))
maxcomp(10, maxd=500, chrs=chrs, nwin=500, winsize=runif(500, 1, 100))
maxcomp(10, maxd=500, chrs=chrs, nwin=600, winsize=runif(600, 1, 100))

maxcomp(10, maxd=500, chrs=chrs, nwin=200, winsize=runif(200, 1, 100))
maxcomp(10, maxd=500, chrs=chrs, nwin=500, winsize=runif(500, 1, 100))
maxcomp(10, maxd=500, chrs=chrs, nwin=600, winsize=runif(600, 1, 100))

maxcomp(10, maxd=200, chrs=chrs, nwin=200, winsize=10)
maxcomp(10, maxd=200, chrs=chrs, nwin=500, winsize=10)
maxcomp(10, maxd=200, chrs=chrs, nwin=600, winsize=10)

maxcomp(10, maxd=200, chrs=chrs, nwin=200, winsize=runif(200, 1, 100))
maxcomp(10, maxd=200, chrs=chrs, nwin=500, winsize=runif(500, 1, 100))
maxcomp(10, maxd=200, chrs=chrs, nwin=600, winsize=runif(600, 1, 100))

###################################################################################################
# Testing the strand-specific nature of clustering.

strcomp <- function(tol=100, maxd=200, ...) {
	stuff <- generateWindows(...)
	strand(stuff) <- sample(c("+", "-", "*"), length(stuff), replace=TRUE)
	combo <- mergeWindows(stuff, tol=tol, max.width=maxd, ignore.strand=FALSE)

	# Checking that each set of strands forms unique IDs.
	out <- split(strand(stuff), combo$id)
	stopifnot(all(sapply(out, FUN=anyDuplicated)==0L))

	# Checking that the strandedness of the output is okay.
	stopifnot(all(strand(combo$region)[combo$id]==strand(stuff)))
	
	# Checking what we get if we set ignore.strand=TRUE.
	combo2 <- mergeWindows(stuff, tol=tol, max.width=maxd, ignore.strand=TRUE)
	stopifnot(all(strand(combo2$region)=="*"))

	# Running separately on each strand, and checking that the boundaries are the same.
	is.forward <- as.logical(strand(stuff)=="+")
	forward <- mergeWindows(stuff[is.forward], tol=tol, max.width=maxd)
	is.reverse <- as.logical(strand(stuff)=="-")
	reverse <- mergeWindows(stuff[is.reverse], tol=tol, max.width=maxd)
	is.unstrand <- as.logical(strand(stuff)=="*")
	unstrand <- mergeWindows(stuff[is.unstrand], tol=tol, max.width=maxd)
	
	strand(forward$region) <- "+"
	strand(reverse$region) <- "-"
	strand(unstrand$region) <- "*"
	if (!identical(c(forward$region, reverse$region, unstrand$region), combo$region)) { 
		stop("mismatch in regions after stranded merging") }

	final.out <- integer(length(stuff))
	final.out[is.forward] <- forward$id
	final.out[is.reverse] <- reverse$id+length(forward$region) 
	final.out[is.unstrand] <- unstrand$id+length(forward$region)+length(reverse$region)
	if (!identical(final.out, combo$id)) { stop("mismatch in IDs after stranded merging") }
	
	return(combo$region)	
}

set.seed(1235213)

strcomp(50, maxd=500, chrs=chrs, nwin=200, winsize=1)
strcomp(50, maxd=500, chrs=chrs, nwin=200, winsize=10)
strcomp(10, maxd=500, chrs=chrs, nwin=200, winsize=100)
	
strcomp(50, chrs=chrs, nwin=500, winsize=1)
strcomp(50, chrs=chrs, nwin=500, winsize=10)
strcomp(10, chrs=chrs, nwin=500, winsize=100)

strcomp(10, maxd=200, chrs=chrs, nwin=200, winsize=runif(200, 1, 100))
strcomp(10, maxd=200, chrs=chrs, nwin=500, winsize=runif(500, 1, 100))
strcomp(10, maxd=200, chrs=chrs, nwin=600, winsize=runif(600, 1, 100))

###################################################################################################
# End.
