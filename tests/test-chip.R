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
		str[current]<-(rbinom(sum(current), 1, 0.5)==1);
	}
	simsam(outfname, names(chromos)[pos.chr], pos.pos, str, chromos);
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
	x<-windowCounts(bamFiles, ext=fraglen, right=right, left=left, spacing=spacing, filter=filter, 
		discard=discard, restrict=restrict)

	# Checking with countOverlaps.
	totals<-integer(length(bamFiles))
	out<-matrix(0L, length(x$region), length(bamFiles))
	for (i in 1:length(bamFiles)) {
        reads <- scanBam(bamFiles[i], param = ScanBamParam(what =c("rname", "strand", "pos", "qwidth")))[[1]]
		read.starts<-ifelse(reads[[2]]=="+", reads[[3]], reads[[3]]+reads[[4]]-fraglen)
		read.ends<-read.starts+fraglen-1L
		frags <- GRanges(reads[[1]], IRanges(read.starts, read.ends))

		if (!is.null(discard)) { frags <- frags[!overlapsAny(GRanges(reads[[1]], IRanges(reads[[3]], reads[[4]]+reads[[3]]-1L)), discard)] }
		if (!is.null(restrict)) { frags <- frags[seqnames(frags) %in% restrict] }
		current<-findOverlaps(x$region, frags)
		out[,i]<-tabulate(queryHits(current), nbins=length(x$region))
		totals[i]<-length(frags)
	}

	if (!identical(out, x$counts)) { stop("mismatch in count matrices") }
	if (!identical(totals, x$totals)) { stop("mismatch in total counts") }

	# Checking the filter.
    x2<-windowCounts(bamFiles, ext=fraglen, right=right, left=left, spacing=spacing, filter=-1, 
			discard=discard, restrict=restrict)
	expected<-expectedRanges(right+left+1L, left, spacing, bamFiles, restrict=restrict)
	if (compare2Ranges(expected, x2$region)) { stop("mismatch in expected and unfiltered regions"); }

	keep<-rowSums(x2$counts)>=filter
	if (!identical(x$counts, x2$count[keep,])) { stop("mismatch in filtered counts"); }
	if (sum(keep)==0 && length(x$region)==0) { } 
	else if (compare2Ranges(x2$region[keep], x$region)) { stop("mismatch in filtered regions"); }

	# Checking MAPQ filtering.
	xx<-windowCounts(bamFiles, minq=200)
	if (!all(xx$totals==0)) { stop("MAPQ filtering failed") }
	return(x$region);
}

# Bin counts.

bincomp <- function(bamFiles, binsize=1000L) {
	binsize<-as.integer(binsize+0.5)
	blah<-windowCounts(bamFiles, bin=binsize)	
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
bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))

# Vanilla comparison.

comp(bamFiles, fraglen=100, spacing=20)
comp(bamFiles, fraglen=200, spacing=50)

bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

# More complex with right arguments.

comp(bamFiles, fraglen=100, right=30, spacing=20)
comp(bamFiles, fraglen=200, left=50, spacing=25)
comp(bamFiles, fraglen=150, right=10, left=10, spacing=30)

# Even more complex, with filtering arguments

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

comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, left=50)

comp(bamFiles, fraglen=200, right=0, spacing=20)
comp(bamFiles, fraglen=200, left=100, spacing=100)

comp(bamFiles, fraglen=100, filter=20)
comp(bamFiles, fraglen=200, filter=40)
comp(bamFiles, fraglen=200, filter=60)

# One more time; sparse across the genome.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))

comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)

bincomp(bamFiles, 1000)
bincomp(bamFiles, 123)
bincomp(bamFiles, 500)

comp(bamFiles, fraglen=100, left=100)
comp(bamFiles, fraglen=200, right=50)

comp(bamFiles, fraglen=200, left=0, spacing=50)
comp(bamFiles, fraglen=200, right=100, spacing=100)

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

comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, right=50)

comp(bamFiles, fraglen=200, right=0, spacing=40)
comp(bamFiles, fraglen=200, right=100, spacing=100)

comp(bamFiles, fraglen=200, filter=100)
comp(bamFiles, fraglen=200, filter=150)

###################################################################################################
# We try to do some crazy spacing.

chromos<-c(chrA=8000, chrB=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))

comp(bamFiles, spacing=max(chromos))

# We also try to do some crazy extension. We should see 1 combination per chromosome, consisting of all reads.

comp(bamFiles, right=max(chromos), left=max(chromos))

###################################################################################################
# Restricted and/or discarded.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))

makeDiscard <- function(ndisc, sizeof) {
	chosen <- sample(length(chromos), ndisc, replace=T)
	chosen.pos <- runif(ndisc, 1, chromos[chosen]-sizeof)
	GRanges(names(chromos)[chosen], IRanges(chosen.pos, chosen.pos+sizeof))
}

comp(bamFiles, fraglen=100, discard=makeDiscard(10, 200))
comp(bamFiles, fraglen=200, discard=makeDiscard(20, 100), restrict="chrA")

comp(bamFiles, fraglen=100, left=100, discard=makeDiscard(50, 20))
comp(bamFiles, fraglen=200, right=50, discard=makeDiscard(10, 200), restrict=c("chrA", "chrB"))

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

