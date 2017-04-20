############################################################
# This tests the profileSites command, to ensure that it's actually giving proper advice.

source("simsam.R")

sdir<-"profile-test"
dir.create(sdir)
outfname <- file.path(sdir, "out")

suppressWarnings(suppressPackageStartupMessages(require(csaw)))

comp <- function(nreads, chromos, ext=100, width=200, res=50, minq=NA, dedup=FALSE, 
		use.strand=FALSE, match.strand=FALSE, final.mode=NA) { 
	# Simulating first.
	bam <- regen(nreads, chromos, outfname)
	windows <- generateWindows(chrs=chromos, winsize=res, nwin=20)
	nwin <- length(windows)
	if (match.strand) { use.strand <- TRUE }
	if (use.strand) { strand(windows) <- sample(c("+", "-", "*"), nwin, replace=TRUE) }

	# Running profileSites.
	xparam <- readParam(minq=minq, dedup=dedup)
    if (!is.na(final.mode)){ ext <- list(ext, final.mode) }
	if (match.strand) { xparam2 <- reform(xparam, forward=NULL) }
	else { xparam2 <- xparam }
	strand.string <- ifelse(use.strand, ifelse(match.strand, "match", "use"), "ignore")
	all.profiles <- profileSites(bam, windows, ext=ext, range=width, average=FALSE, param=xparam2, strand=strand.string)

    # Checking names.
    if (!identical(colnames(all.profiles), as.character((-width):width))) {
        stop("column names are not as expected")
    }

    # Checking the normalization/averaging options are good.
	none <- profileSites(bam, windows, ext=ext, range=width, param=xparam2, strand=strand.string)
    ref.none <- colMeans(all.profiles)
	if (!isTRUE(all.equal(ref.none, none))) { stop("unnormalized coverage profiles don't match up with manual average") }
	total <- profileSites(bam, windows, ext=ext, range=width, param=xparam2, strand=strand.string, normalize="total")
    ref.total <- colMeans(all.profiles/rowSums(all.profiles))
	if (!isTRUE(all.equal(ref.total, total))) { stop("total-normalized coverage profiles don't match up with manual average") }
	max <- profileSites(bam, windows, ext=ext, range=width, param=xparam2, strand=strand.string, normalize="max")
    ref.max <- colMeans(all.profiles/apply(all.profiles, 1, max))
	if (!isTRUE(all.equal(ref.max, max))) { stop("max-normalized coverage profiles don't match up with manual average") }

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
	for (x in seq_len(nwin)) {
		curchr <- as.character(seqnames(windows[x]))
		relevant <- relevant.start[x]:relevant.end[x]
		valid <- relevant > 0L & relevant <= chromos[[curchr]]
				
		# Using reverse coverage if match.strand is TRUE.
		if (match.strand && reverse[x]) { 
			chosen <- totally.reverse			
		} else {
			chosen <- totally
		}
        cur.prof <- as.integer(chosen[[curchr]][relevant[valid]])
        if (!identical(cur.prof, unname(all.profiles[x,valid]))) { 
            stop("region-specific coverage profile doesn't match up") 
        } else if (!all(all.profiles[x,!valid]==0)) {
            stop("non-zero entries for invalid coverage values")
        }
	}

	# Evaluating the two methods.
	return(head(none))
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

comp(nreads, chromos, res=20, dedup=TRUE)
comp(nreads, chromos, res=100, minq=100)
comp(nreads, chromos, width=100, dedup=TRUE)
comp(nreads, chromos, width=500, minq=100)

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

comp(nreads, chromos, ext=50, final.mode=50)
comp(nreads, chromos, ext=50, final.mode=100)
comp(nreads, chromos, ext=50, final.mode=20)

############################################################
# Cleaning up.

unlink(sdir, recursive=TRUE)

############################################################
