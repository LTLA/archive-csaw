###################################################################################################
# This script tests the ChIP capabilities of the 'csaw' package.

suppressWarnings(suppressPackageStartupMessages(library(csaw)))
source("simsam.R")

expectedRanges <- function(width, offset, spacing, bam.files, restrict=NULL) {
	spacing<-as.integer(spacing)
	width<-as.integer(width)
	offset<-as.integer(offset)
	chrs<-scanBamHeader(bam.files[1])[[1]][[1]]
	output<-list()
	if (!is.null(restrict)) { chrs <- chrs[names(chrs) %in% restrict] }

	for (x in names(chrs)) {
		multiples<-ceiling(chrs[[x]]/spacing)
		all.starts<-0:multiples*spacing+1L-offset
		all.ends<-all.starts+width-1L
		all.starts<-pmax(all.starts, 1L)
		all.ends<-pmin(all.ends, chrs[[x]])

		keep <- all.starts <= chrs[[x]] & all.ends > 0L
		gr <-GRanges(x, IRanges(all.starts[keep], all.ends[keep]))
		keep <- !GenomicRanges::duplicated(gr) 
		output[[length(output)+1L]]<- gr[keep]
	}
	output<-suppressWarnings(do.call(c, output))
	return(output)
}

compare2Ranges <- function(left, right) {
	left<-sort(left)
	right<-sort(right)
	if (length(left)!=length(right)) { return(TRUE) }
	if (any(as.character(seqnames(left))!=as.character(seqnames(right))) ||
	    any(start(left)!=start(right)) || any((end(left)!=end(right)))) { return(TRUE); }
	return(FALSE)
}

# We also set up a comparison function between the windowCount function and its countOverlaps equivalent.

comp <- function(bamFiles, fraglen=200, right=0, left=0, spacing=20, filter=-1, discard=GRanges(), restrict=NULL, forward=NA) {
	if (length(fraglen)==1L) { 
		fraglen <- rep(fraglen, length.out=length(bamFiles))
		remainder <- integer(length(bamFiles))
	} else {
		final.out <- as.integer(mean(fraglen))
		remainder <- as.integer((final.out - fraglen)/2)
		fraglen <- makeExtVector(fraglen, final.out)
	}
	chrlens <- csaw:::.activeChrs(bamFiles, NULL)
	
	for (type in 1:3) {
		if (type==1) {
			dedup<- FALSE
			minq <- 0
		} else if (type==2) {
			dedup <- TRUE
			minq <- 0
		} else if (type==3) {
			dedup <- FALSE
			minq <- 100
		}
		x<-windowCounts(bamFiles, ext=fraglen, width=right+left+1, shift=left, spacing=spacing, filter=filter, 
			param=readParam(discard=discard, restrict=restrict, minq=minq, dedup=dedup, forward=forward))

		# Checking with countOverlaps.
		totals<-integer(length(bamFiles))
		out<-matrix(0L, nrow(x), length(bamFiles))
		for (i in 1:length(bamFiles)) {
	        reads <- scanBam(bamFiles[i], param = ScanBamParam(what =c("rname", "strand", "pos", "qwidth", "mapq"), 
						flag=scanBamFlag(isDuplicate=ifelse(dedup, FALSE, NA), 
							isMinusStrand=ifelse(is.na(forward), NA, !forward))))[[1]]
			keep <- reads$mapq >= minq
			reads$mapq <- NULL
			if (any(keep)) {
				reads$pos <- reads$pos[keep]
				reads$strand <- reads$strand[keep]
				reads$rname <- reads$rname[keep]
				reads$qwidth <- reads$qwidth[keep]
			}

			read.starts<-ifelse(reads$strand=="+", reads$pos, reads$pos+reads$qwidth-fraglen[i])
			read.ends<-read.starts+fraglen[i]-1L
			read.starts <- read.starts - remainder[i]
			read.ends <- read.ends + remainder[i]
			if (length(read.starts)) { read.starts <- pmin(read.starts, chrlens[reads$rname]) }
			if (length(read.ends)) { read.ends <- pmax(1L, read.ends) }
			frags <- GRanges(reads$rname, IRanges(read.starts, read.ends))

			# Discarding. No variable read lengths here, so no need to use alignment width.			
			if (length(discard)) { frags <- frags[!overlapsAny(GRanges(reads$rname, IRanges(reads$pos, reads$pos+reads$qwidth-1L)), discard, type="within")] }
			if (!is.null(restrict)) { frags <- frags[seqnames(frags) %in% restrict] }
			out[,i]<-countOverlaps(rowRanges(x), frags)
			totals[i]<-length(frags)
		}
		
		curcounts <- assay(x)	
		attributes(curcounts)$dimnames <- NULL
		if (!identical(out, curcounts)) { stop("mismatch in count matrices") }
		if (!identical(totals, x$totals)) { stop("mismatch in total counts") }

		# Checking the filter. We need to do this separately as the check above is not filter-aware.
		if (filter==-1) {
			x2 <- x
		} else {
	    	x2<-windowCounts(bamFiles, ext=fraglen, width=right+left+1, shift=left, spacing=spacing, filter=-1, 
				param=readParam(discard=discard, restrict=restrict, dedup=dedup, minq=minq, forward=forward))
			keep<-rowSums(assay(x2))>=filter
			if (!identical(assay(x), assay(x2)[keep,])) { stop("mismatch in filtered counts") }
			if (sum(keep)==0 && nrow(x)==0) { } 
			else if (compare2Ranges(rowRanges(x2)[keep], rowRanges(x))) { stop("mismatch in filtered regions") }
		}
		if (type==1) { 
			expected<-expectedRanges(right+left+1L, left, spacing, bamFiles, restrict=restrict)
			if (compare2Ranges(expected, rowRanges(x2))) { stop("mismatch in expected and unfiltered regions") }
		}
	}

	return(rowRanges(x))
}

# Bin count checker.

bincomp <- function(bamFiles, binsize=1000L) {
	binsize<-as.integer(binsize+0.5)
	blah<-windowCounts(bamFiles, width=binsize, bin=TRUE)	
	if (!identical(blah$totals, as.integer(colSums(assay(blah)) ) ) ) { stop("totals do not match up") }
	expected<-expectedRanges(binsize, 0L, binsize, bamFiles)
	chrs<-scanBamHeader(bamFiles)[[1]][[1]]

	# Counting reads into bins.
	total.out<-list()
	regions <- rowRanges(blah)
    for (x in runValue(seqnames(regions))) {
    	out<-matrix(0L, ceiling(chrs[[x]]/binsize), length(bamFiles))
		for (i in 1:length(bamFiles)) { 
			reads <- scanBam(bamFiles[i], param = ScanBamParam(what =c("rname", "strand", "pos", "qwidth"),
						which=GRanges(x, IRanges(1, chrs[[x]]))))[[1]]
			pos<-ifelse(reads[[2]]=="+", reads[[3]], reads[[3]]+reads[[4]]-1L)
			pos<-pmin(pos, chrs[[x]])
			out[,i]<-tabulate((pos-1L)/binsize+1L, nbins=nrow(out))
		}

		keep<-rowSums(out)>0L
		testing <- regions[seqnames(regions)==x]
		expect<-expected[seqnames(expected)==x]
		if (compare2Ranges(testing, expect[keep])) { 
			stop("bins do not match up after filtering") }
		total.out[[x]]<-out[keep,]
   	}

	# Comparing counts.
	total.out<-do.call(rbind, total.out)
	curcounts <- assay(blah)
	attributes(curcounts)$dimnames <- NULL
	if (!identical(total.out, curcounts)) { stop("binned counts do not match up") }
	return(blah$totals)
}

###################################################################################################
# Setting up some variables to do the comparison.

dir<-"chip-test";
dir.create(dir);

set.seed(2123)
chromos<-c(chrA=10000, chrB=5000)

# Vanilla comparison.

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, spacing=20)
comp(bamFiles, fraglen=200, spacing=50)
comp(bamFiles, fraglen=c(100, 200), spacing=50)
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

# More complex with right arguments.

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, right=30, spacing=20)
comp(bamFiles, fraglen=200, left=20, spacing=25)
comp(bamFiles, fraglen=c(100, 200), left=20, spacing=25)
comp(bamFiles, fraglen=150, right=10, left=10, spacing=30)
comp(bamFiles, fraglen=150, right=-5, left=15)

# Even more complex, with filtering arguments

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=5)
comp(bamFiles, fraglen=100, filter=10)
comp(bamFiles, fraglen=200, filter=15)
comp(bamFiles, fraglen=c(100, 200), filter=15)

# And again, with a different chromosome set-up.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)
comp(bamFiles, fraglen=c(100,200))
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, left=10)
comp(bamFiles, fraglen=200, right=-30, left=40, spacing=50)
comp(bamFiles, fraglen=200, left=70, spacing=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=20)
comp(bamFiles, fraglen=200, filter=40)
comp(bamFiles, fraglen=200, filter=60)

# One more time; sparse across the genome, but three files.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)
comp(bamFiles, fraglen=c(100, 150, 200))
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, left=10)
comp(bamFiles, fraglen=200, right=50)
comp(bamFiles, fraglen=200, right=-30, left=40, spacing=50)
comp(bamFiles, fraglen=200, right=70, spacing=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=0)
comp(bamFiles, fraglen=200, filter=1)
comp(bamFiles, fraglen=200, right=50, filter=2)

###################################################################################################
# We test with four files. Oh, the humanity.

bamFiles<-sapply(1:4, FUN=function(i) { regen(3000, chromos, file.path(dir, paste("A", i, sep=""))) })
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)
comp(bamFiles, fraglen=c(100, 120, 140, 160))
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

bamFiles<-sapply(1:4, FUN=function(i) { regen(3000, chromos, file.path(dir, paste("A", i, sep=""))) })
comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, right=50)
comp(bamFiles, fraglen=200, left=20, right=-15, spacing=40)
comp(bamFiles, fraglen=200, left=50, right=100, spacing=100)
comp(bamFiles, fraglen=200, filter=100)
comp(bamFiles, fraglen=200, filter=150)

###################################################################################################
# We try to do some crazy spacing and some crazy extension. We should see 1 combination per chromosome, consisting of all reads.

chromos<-c(chrA=8000, chrB=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, spacing=max(chromos))
comp(bamFiles, right=max(chromos))
comp(bamFiles, right=max(chromos), spacing=max(chromos))

###################################################################################################
# Restricted and/or discarded and/or strand-specific.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, discard=makeDiscard(10, 200, chromos))
comp(bamFiles, fraglen=200, discard=makeDiscard(20, 100, chromos), restrict="chrA")
comp(bamFiles, fraglen=100, left=30, spacing=50, discard=makeDiscard(50, 20, chromos))
comp(bamFiles, fraglen=200, right=50, discard=makeDiscard(10, 200, chromos), restrict=c("chrA", "chrB"))

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=200, left=20, spacing=50, discard=makeDiscard(20, 200, chromos))
comp(bamFiles, fraglen=200, right=100, spacing=100, discard=makeDiscard(10, 1000, chromos))
comp(bamFiles, fraglen=100, filter=0, discard=makeDiscard(10, 500, chromos))
comp(bamFiles, fraglen=200, filter=1, discard=makeDiscard(5, 1000, chromos), restrict=c("chrC", "chrA"))
comp(bamFiles, fraglen=200, right=50, filter=2, discard=makeDiscard(20, 100, chromos))

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, forward=TRUE)
comp(bamFiles, fraglen=100, forward=FALSE, discard=makeDiscard(10, 200, chromos))
comp(bamFiles, fraglen=200, left=20, spacing=50, forward=TRUE)
comp(bamFiles, fraglen=200, left=20, spacing=50, forward=FALSE, discard=makeDiscard(20, 200, chromos))

###################################################################################################
# Cleaning up.

unlink(dir, recursive=TRUE)

###################################################################################################
# End.

