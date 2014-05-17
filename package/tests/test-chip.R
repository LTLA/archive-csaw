###################################################################################################
# This script tests the ChIP capabilities of the 'csaw' package.

suppressPackageStartupMessages(library(csaw))
source("simsam.R")

# First, we set up some functions to generate some random SAM files.

regen <- function(nreads, chromos, outfname) {
	pos.chr<-sample(length(chromos), nreads, replace=TRUE)
	pos.pos<-rep(0, nreads)
	str<-rep(0, nreads)
	for (i in 1:length(chromos)) {
		current<-pos.chr==i
		pos.pos[current]<-round(runif(sum(current), 1, chromos[i]))
		str[current]<-(rbinom(sum(current), 1, 0.5)==1)
	}
	isdup <- rbinom(nreads, 1, 0.8)==0L
	mapq <- round(runif(nreads, 50, 199))
	simsam(outfname, names(chromos)[pos.chr], pos.pos, str, chromos, is.dup=isdup, mapq=mapq)
}

expectedRanges <- function(width, offset, spacing, bam.files, restrict=NULL) {
	spacing<-as.integer(spacing)
	width<-as.integer(width)
	offset<-as.integer(offset)
	chrs<-scanBamHeader(bam.files[1])[[1]][[1]]
	output<-list()
	if (!is.null(restrict)) { chrs <- chrs[names(chrs) %in% restrict] }

	for (x in names(chrs)) {
		multiples<-as.integer((chrs[[x]]-1L)/spacing)
		all.starts<-0:multiples*spacing+1L-offset
		all.ends<-all.starts+width-1L
		all.starts<-pmax(all.starts, 1L)
		all.ends<-pmin(all.ends, chrs[[x]])

		gr <-GRanges(x, IRanges(all.starts, all.ends))
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

comp <- function(bamFiles, fraglen=200, right=0, left=0, spacing=20, filter=-1, discard=NULL, restrict=NULL) {
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
			discard=discard, restrict=restrict, minq=minq, dedup=dedup)

		# Checking with countOverlaps.
		totals<-integer(length(bamFiles))
		out<-matrix(0L, length(x$region), length(bamFiles))
		for (i in 1:length(bamFiles)) {
	        reads <- scanBam(bamFiles[i], param = ScanBamParam(what =c("rname", "strand", "pos", "qwidth", "mapq"), 
						flag=scanBamFlag(isDuplicate=ifelse(dedup, FALSE, NA))))[[1]]
			keep <- reads$mapq >= minq
			reads$mapq <- NULL
			if (any(keep)) {
				reads$pos <- reads$pos[keep]
				reads$strand <- reads$strand[keep]
				reads$rname <- reads$rname[keep]
				reads$qwidth <- reads$qwidth[keep]
			}

			read.starts<-ifelse(reads$strand=="+", reads$pos, reads$pos+reads$qwidth-fraglen)
			read.ends<-read.starts+fraglen-1L
			frags <- GRanges(reads$rname, IRanges(read.starts, read.ends))

			# Discarding. No variable read lengths here, so no need to use alignment width.			
			if (!is.null(discard)) { frags <- frags[!overlapsAny(GRanges(reads$rname, IRanges(reads$pos, reads$pos+reads$qwidth-1L)), discard, type="within")] }
			if (!is.null(restrict)) { frags <- frags[seqnames(frags) %in% restrict] }
			out[,i]<-countOverlaps(x$region, frags)
			totals[i]<-length(frags)
		}

		if (!identical(out, x$counts)) { stop("mismatch in count matrices") }
		if (!identical(totals, x$totals)) { stop("mismatch in total counts") }

		# Checking the filter. We need to do this separately as the check above is not filter-aware.
		if (filter==-1) {
			x2 <- x
		} else {
	    	x2<-windowCounts(bamFiles, ext=fraglen, width=right+left+1, shift=left, spacing=spacing, filter=-1, 
				discard=discard, restrict=restrict, dedup=dedup, minq=minq)
			keep<-rowSums(x2$counts)>=filter
			if (!identical(x$counts, x2$count[keep,])) { stop("mismatch in filtered counts") }
			if (sum(keep)==0 && length(x$region)==0) { } 
			else if (compare2Ranges(x2$region[keep], x$region)) { stop("mismatch in filtered regions") }
		}
		if (type==1) { 
			expected<-expectedRanges(right+left+1L, left, spacing, bamFiles, restrict=restrict)
			if (compare2Ranges(expected, x2$region)) { stop("mismatch in expected and unfiltered regions") }
		}
	}

	return(x$region);
}

# Bin count checker.

bincomp <- function(bamFiles, binsize=1000L) {
	binsize<-as.integer(binsize+0.5)
	blah<-windowCounts(bamFiles, width=binsize, bin=TRUE)	
	if (!identical(blah$totals, as.integer(colSums(blah$counts)+0.5))) { stop("totals do not match up") }
	expected<-expectedRanges(binsize, 0L, binsize, bamFiles)
	chrs<-scanBamHeader(bamFiles)[[1]][[1]]

	# Counting reads into bins.
	total.out<-list()
    for (x in runValue(seqnames(blah$region))){
    	out<-matrix(0L, ceiling(chrs[[x]]/binsize), length(bamFiles))
		for (i in 1:length(bamFiles)) { 
			reads <- scanBam(bamFiles[i], param = ScanBamParam(what =c("rname", "strand", "pos", "qwidth"),
						which=GRanges(x, IRanges(1, chrs[[x]]))))[[1]]
			pos<-ifelse(reads[[2]]=="+", reads[[3]], reads[[3]]+reads[[4]]-1L)
			pos<-pmin(pos, chrs[[x]])
			out[,i]<-tabulate((pos-1L)/binsize+1L, nbins=nrow(out))
		}

		keep<-rowSums(out)>0L
		testing<-blah$region[seqnames(blah$region)==x]
		expect<-expected[seqnames(expected)==x]
		if (compare2Ranges(testing, expect[keep])) { 
			stop("bins do not match up after filtering") }
		total.out[[x]]<-out[keep,]
   	}
	total.out<-do.call(rbind, total.out)
	if (!identical(total.out, blah$counts)) { stop("binned counts do not match up") }
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
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

# More complex with right arguments.

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, right=30, spacing=20)
comp(bamFiles, fraglen=200, left=50, spacing=25)
comp(bamFiles, fraglen=150, right=10, left=10, spacing=30)

# Even more complex, with filtering arguments

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=5)
comp(bamFiles, fraglen=100, filter=10)
comp(bamFiles, fraglen=200, filter=15)

# And again, with a different chromosome set-up.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, left=50)
comp(bamFiles, fraglen=200, right=0, spacing=20)
comp(bamFiles, fraglen=200, left=100, spacing=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=20)
comp(bamFiles, fraglen=200, filter=40)
comp(bamFiles, fraglen=200, filter=60)

# One more time; sparse across the genome, but three files.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, left=100)
comp(bamFiles, fraglen=200, right=50)
comp(bamFiles, fraglen=200, left=0, spacing=50)
comp(bamFiles, fraglen=200, right=100, spacing=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=0)
comp(bamFiles, fraglen=200, filter=1)
comp(bamFiles, fraglen=200, right=50, filter=2)

###################################################################################################
# We test with four files. Oh, the humanity.

bamFiles<-sapply(1:4, FUN=function(i) { regen(3000, chromos, file.path(dir, paste("A", i, sep=""))) })
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)
bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

bamFiles<-sapply(1:4, FUN=function(i) { regen(3000, chromos, file.path(dir, paste("A", i, sep=""))) })
comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, right=50)
comp(bamFiles, fraglen=200, right=0, spacing=40)
comp(bamFiles, fraglen=200, right=100, spacing=100)
comp(bamFiles, fraglen=200, filter=100)
comp(bamFiles, fraglen=200, filter=150)

###################################################################################################
# We try to do some crazy spacing and some crazy extension. We should see 1 combination per chromosome, consisting of all reads.

chromos<-c(chrA=8000, chrB=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, spacing=max(chromos))
comp(bamFiles, right=max(chromos), left=max(chromos))

###################################################################################################
# Restricted and/or discarded.

makeDiscard <- function(ndisc, sizeof) {
	chosen <- sample(length(chromos), ndisc, replace=T)
	chosen.pos <- runif(ndisc, 1, chromos[chosen]-sizeof)
	reduce(GRanges(names(chromos)[chosen], IRanges(chosen.pos, chosen.pos+sizeof)))
}
chromos<-c(chrA=5000, chrB=5000, chrC=8000)

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, discard=makeDiscard(10, 200))
comp(bamFiles, fraglen=200, discard=makeDiscard(20, 100), restrict="chrA")
comp(bamFiles, fraglen=100, left=100, discard=makeDiscard(50, 20))
comp(bamFiles, fraglen=200, right=50, discard=makeDiscard(10, 200), restrict=c("chrA", "chrB"))

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=200, left=0, spacing=50, discard=makeDiscard(20, 200))
comp(bamFiles, fraglen=200, right=100, spacing=100, discard=makeDiscard(10, 1000))
comp(bamFiles, fraglen=100, filter=0, discard=makeDiscard(10, 500))
comp(bamFiles, fraglen=200, filter=1, discard=makeDiscard(5, 1000), restrict=c("chrC", "chrA"))
comp(bamFiles, fraglen=200, right=50, filter=2, discard=makeDiscard(20, 100))

###################################################################################################
# Cleaning up.

unlink(dir, recursive=TRUE)

###################################################################################################
# End.

