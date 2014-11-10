############################################################
# This tests the profileSites command, to ensure that it's actually giving proper advice.

source("simsam.R")

sdir<-"profile-test"
dir.create(sdir)
outfname <- file.path(sdir, "out")

suppressWarnings(suppressPackageStartupMessages(require(csaw)))

comp <- function(nreads, chromos, ext=100, width=200, res=50, weight=TRUE, minq=NA, dedup=FALSE) { 
	# Simulating first.
	bam <- regen(nreads, chromos, outfname)
	win.data <- generateWindows(chrs=chromos, winsize=res, nwin=20)
	windows <- win.data$region
	nwin <- length(windows)

	# Running profileSites.
	xparam <- readParam(minq=minq, dedup=dedup)
	if (weight) {
		by.win <- regionCounts(bam, windows, ext=ext, param=xparam)
		metric <- rowSums(assay(by.win))
	} else {
		metric <- rep(1, nwin)
	}
	observed <- profileSites(bam, windows, ext=ext, range=width, param=xparam, weight=metric)

	# Running the reference analysis.
	totally <- list()
	for (chr in names(chromos)) {
		out <- extractReads(GRanges(chr, IRanges(1, chromos[[chr]])), bam, param=xparam)
		totally[[chr]] <- suppressWarnings(resize(out, width=ext))
	} 
	totes.prof <- 0
	for (x in 1:nwin) {
		curwin <- windows[x]
		all.reads <- totally[[as.character(seqnames(curwin))]]
		dist.back <- start(all.reads) - end(curwin)
		dist.front <- start(curwin) - end(all.reads)
		totes.prof <- totes.prof + (tabulate(dist.back, nbin=width) + tabulate(dist.front, nbin=width))/metric[x]
	}

	# Evaluating the two methods.
	reference <- totes.prof/nwin/2
	if (length(reference)!=length(observed)) { stop("vectors are of differing lengths") }
	if (any(abs(reference - observed) > (reference+1e-3)*1e-6)) { stop("coverage profiles don't match up") }
	return(head(observed))
}

############################################################
# Fairly hefty simulations are necessary here.

set.seed(123123)
nreads <- 5000
chromos <- c(chrA=10000, chrB=5000)
comp(nreads, chromos)
comp(nreads, chromos, minq=100)
comp(nreads, chromos, dedup=TRUE)

comp(nreads, chromos, ext=50)
comp(nreads, chromos, ext=200)
comp(nreads, chromos, width=100)
comp(nreads, chromos, width=500)

comp(nreads, chromos, res=20)
comp(nreads, chromos, res=20, width=100)
comp(nreads, chromos, res=100)
comp(nreads, chromos, res=100, width=500)

comp(nreads, chromos, res=20, weight=FALSE)
comp(nreads, chromos, res=100, weight=FALSE)
comp(nreads, chromos, weight=FALSE)
comp(nreads, chromos, weight=FALSE)

############################################################
# Cleaning up.

unlink(sdir, recursive=TRUE)

############################################################
