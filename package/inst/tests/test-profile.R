############################################################
# This tests the profileSites command, to ensure that it's actually giving proper advice.

source("simsam.R")

sdir<-"profile-test"
dir.create(sdir)
outfname <- file.path(sdir, "out")

suppressWarnings(suppressPackageStartupMessages(require(csaw)))

comp <- function(nreads, chromos, ext=100, width=200, res=50, weight=TRUE, minq=NA, dedup=FALSE, 
		use.strand=FALSE, match.strand=FALSE, final.mode=NA) { 
	# Simulating first.
	bam <- regen(nreads, chromos, outfname)
	windows <- generateWindows(chrs=chromos, winsize=res, nwin=20)
	nwin <- length(windows)
	if (match.strand) { use.strand <- TRUE }
	if (use.strand) { strand(windows) <- sample(c("+", "-", "*"), nwin, replace=TRUE) }

	# Running profileSites.
	xparam <- readParam(minq=minq, dedup=dedup)
	if (weight) {
		suppressWarnings(by.win <- regionCounts(bam, windows, ext=ext, param=xparam))
		metric <- rowSums(assay(by.win))
	} else {
		metric <- rep(1, nwin)
	}
	ext <- makeExtVector(ext, final.mode)
	
	if (match.strand) { xparam2 <- reform(xparam, forward=NULL) }
	else { xparam2 <- xparam }
	strand.string <- ifelse(use.strand, ifelse(match.strand, "match", "use"), "ignore")
	observed <- profileSites(bam, windows, ext=ext, range=width, param=xparam2, strand=strand.string, weight=1/metric)

	# Checking it's the same if we do a weighted average.
	all.profiles <- profileSites(bam, windows, ext=ext, range=width, average=FALSE, param=xparam2, strand=strand.string)
	other.observed <- colMeans(all.profiles/metric)
	if (length(other.observed)!=length(observed)) { stop("vectors are of differing lengths against manual average") }
	if (any(abs(other.observed - observed) > (other.observed+1e-3)*1e-6)) { stop("coverage profiles don't match up with manual average") }

	# Running the reference analysis.
	totally <- totally.reverse <- list()
	for (chr in names(chromos)) {
		out <- extractReads(bam, GRanges(chr, IRanges(1, chromos[[chr]])), param=xparam, ext=ext)
		if (!match.strand) { 
			totally[[chr]] <- coverage(ranges(out), width=chromos[[chr]]) 
		} else {
			rev.read <- strand(out)=="-"
			totally[[chr]] <- coverage(ranges(out)[!rev.read], width=chromos[[chr]]) 
			totally.reverse[[chr]] <- coverage(ranges(out)[rev.read], width=chromos[[chr]]) 
		}
	} 

	relevant.start <- start(windows) - width
	relevant.end <- start(windows) + width
	if (use.strand) { 
		reverse <- as.logical(strand(windows)=="-")
		relevant.start[reverse] <- end(windows[reverse]) + width # Automatic reversal for reverse-stranded regions.
		relevant.end[reverse] <- end(windows[reverse]) - width
	}
	totes.prof <- integer(width*2+1)
	for (x in 1:nwin) {
		curchr <- as.character(seqnames(windows[x]))
		relevant <- relevant.start[x]:relevant.end[x]
		valid <- relevant > 0L & relevant <= chromos[[curchr]]
				
		# Using reverse coverage if match.strand is TRUE.
		if (match.strand && reverse[x]) { 
			chosen <- totally.reverse			
		} else {
			chosen <- totally
		}
		totes.prof[valid] <- totes.prof[valid] + as.integer(chosen[[curchr]][relevant[valid]])/metric[x]
	}

	# Evaluating the two methods.
	reference <- totes.prof/nwin
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

# Checking effect of stranded options.

comp(nreads, chromos, res=20, use.strand=TRUE)
comp(nreads, chromos, res=100, use.strand=TRUE)
comp(nreads, chromos, use.strand=TRUE)
comp(nreads, chromos, use.strand=TRUE)

comp(nreads, chromos, res=20, match.strand=TRUE)
comp(nreads, chromos, res=100, match.strand=TRUE)
comp(nreads, chromos, match.strand=TRUE)
comp(nreads, chromos, match.strand=TRUE)

# Just exercising the multi-fragment length options here.
comp(nreads, chromos, ext=50, final.mode=NULL)
comp(nreads, chromos, ext=50, final.mode=100)
comp(nreads, chromos, ext=50, final.mode=20)

############################################################
# Cleaning up.

unlink(sdir, recursive=TRUE)

############################################################
