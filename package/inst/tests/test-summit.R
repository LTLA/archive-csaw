############################################################
# This tests the profileSummit command, to ensure that it's actually giving proper advice.

source("simsam.R")

sdir<-"summit-test"
dir.create(sdir)
outfname <- file.path(sdir, "out")

suppressWarnings(suppressPackageStartupMessages(require(csaw)))

comp <- function(nreads, chromos, ext=100, width=200, res=50, min.depth=20, minq=NA, dedup=FALSE) { 
	# Simulating first.
	bam <- regen(nreads, chromos, outfname)

	# Running the profileSummits.
	xparam <- readParam(minq=minq, dedup=dedup)
	out <- profileSummit(bam, ext=ext, width=width, res=res, min.depth=min.depth/res, param=xparam)

	# Running the reference analysis.
	res <- as.integer(res)
	actual.width <- as.integer(width/res)
	min.depth <- as.integer(min.depth)
	profile <- num <- numeric(actual.width)

	data <- windowCounts(bam, param=xparam, ext=ext, width=res, spacing=res)
	for (x in names(chromos)) {
		nbins <- ceiling(chromos[[x]]/res)
		track <- integer(nbins)
		current <- as.logical(seqnames(rowData(data))==x)
		track[(start(rowData(data))[current] - 1L)/res + 1L] <- assay(data)[current,]

		# Identifying the local maxima.
		all.max <- list()
		for (y in 1:length(track)) {
			if (track[y]<min.depth) { next }

			lower <- max(1L, y-actual.width)
			upper <- min(y+actual.width, length(track))
			cur.track <- track[lower:upper]

			if (sum(cur.track>=track[y])==1L) {
				all.max[[length(all.max)+1L]] <- y
			}
		}

		# Expanding windows and evaluating.
		for (y in all.max) {
			new.core <- y*res
			new.starts <- new.core - 0:actual.width * res + 1L - res
			no.start <- new.starts <= 0L
			new.ends <- new.core + 0:actual.width*res
			no.end <- new.ends > chromos[[x]]

			new.ends[no.end] <- chromos[[x]]
			new.starts[no.start] <- 1L
			new.ranges <- GRanges(x, IRanges(new.starts, new.ends))
			datax <- regionCounts(bam, new.ranges, ext=ext, param=xparam)

			profile <- profile + diff(assay(datax)[,1]) / track[y]
			num <- num + as.integer(!no.start[-1]) + as.integer(!no.end[-1])
		}
	}

	# Evaluating the two methods.
	reference <- profile/num * 2
	if (length(reference)!=length(out$coverage)) { stop("vectors are of differing lengths") }
	if (!identical(is.na(reference), is.na(out$coverage))) { stop("NA values are not identical") }
	if (any(abs(reference - out$coverage) > (reference+1e-3)*1e-6, na.rm=TRUE)) { stop("summit profiles don't match up") }
	return(head(out$coverage))
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

comp(nreads, chromos, res=20, min.depth=10)
comp(nreads, chromos, res=100, min.depth=50)
comp(nreads, chromos, min.depth=10)
comp(nreads, chromos, min.depth=50)

############################################################
# Cleaning up.

unlink(sdir, recursive=TRUE)

############################################################
