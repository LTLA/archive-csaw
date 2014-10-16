###################################################################################################
# We test the correlateChIP function in 'csaw' against an equivalent version in R.

source("simsam.R")

fdir<-"ccf-test"
dir.create(fdir)
outfname<-file.path(fdir, "out")

suppressWarnings(suppressPackageStartupMessages(library(csaw)))

manualcor<-function(bamx, n, cross, minq=0, dedup=FALSE) { 
	chromos<-scanBamHeader(bamx)[[1]][[1]]
	out<-0
	total<-0
	for (chr in names(chromos)) {
		clen<-chromos[[chr]]
		param <- ScanBamParam(what=c("pos", "qwidth", "strand", "mapq"), which=GRanges(chr, IRanges(1, clen)),
			flag=scanBamFlag(isDuplicate=ifelse(dedup, FALSE, NA)))
		reads<-list()
		for (b in bamx) {
			new.reads <- scanBam(b, param = param)[[1]]	
			keep <- new.reads$mapq >= minq
			reads$str<-c(reads$str, new.reads$str[keep]=="+")
			reads$qwidth<-c(reads$qwidth, new.reads$qwidth[keep])
			reads$pos<-c(reads$pos, new.reads$pos[keep])
		}
		f<-r<-rep(0, clen)

		fx<-table(reads$pos[reads$str])
		f[as.integer(names(fx))]<-as.integer(fx)

		rx<-table(pmin(clen, reads$pos[!reads$str]+reads$qwidth[!reads$str]))
		r[as.integer(names(rx))]<-as.integer(rx)

		# Autocorrelations, if not cross-correaltions, so we just fuse them together.
		if (!cross) {
			f<-f+r
			r<-f
		}

		nreads<-length(reads$pos)
		out<-out+nreads*sapply(0:n, FUN=function(i){ 
			if (i>=length(f)-1L || i>=length(r)-1L) { return(0) }
			fr<-f[1:(length(f)-i)]
			rr<-r[(i+1):length(r)]
			if (sd(fr)==0 || sd(rr)==0) { return(0) }
			cor(fr, rr)
		})
		total<-total+nreads
	}
	out/total
}

comp<-function(bamFiles, n, cross=TRUE) {
	precision<-1e-8
    for (type in 1:3) {
        if (type==1) {
            dedup <- FALSE
            minq <- 100
        } else if (type==2) {
            dedup <- TRUE
            minq <- 0
        } else if (type==3) {
            dedup<- FALSE
            minq <- 0
        }
		out<-manualcor(bamFiles, n, cross=cross, minq=minq, dedup=dedup)
		out2<-correlateReads(bamFiles, n, cross=cross, readParam(minq=minq, dedup=dedup))
		if (length(out)!=length(out2)) { stop("mismatch in length of output vector"); }
		if (any( abs((out-out2)/(abs(out)+precision)) > precision ))  { stop("mismatch in correlation coefficients"); }
	}
	head(out)
}

###################################################################################################
# Testing with some data.

set.seed(10);
bamFile <- regen(1000, c(chrA=10000), outfname)
comp(bamFile, 50)
comp(bamFile, 100)

# And again...

bamFile <- regen(1000, c(chrA=10000), outfname)
comp(bamFile, 50)
comp(bamFile, 100)

# Repeating with more reads.

bamFile <- regen(2000, c(chrA=10000), outfname)
comp(bamFile, 50)
comp(bamFile, 100)

# Trying it out with multiple chromosomes.

bamFile <- regen(5000, c(chrA=10000, chrB=5000), outfname)
comp(bamFile, 50)
comp(bamFile, 100)

# And again, with more reads.

bamFile <- regen(10000, c(chrA=10000, chrB=5000), outfname)
comp(bamFile, 50)
comp(bamFile, 100)

# Trying it out with multiple BAM files.

bamFiles <- c(regen(500, c(chrA=1000, chrB=500), outfname), regen(500, c(chrA=1000, chrB=500), outfname))
comp(bamFiles, 50)
comp(bamFiles, 100)

# And again, with more reads.

bamFiles <- c(regen(5000, c(chrA=10000, chrB=5000), outfname), regen(5000, c(chrA=10000, chrB=5000), outfname))
comp(bamFiles, 50)
comp(bamFiles, 100)

###################################################################################################
# Repeating; but this time, looking at autocorrelations.

set.seed(1034785)
bamFile <- regen(1000, c(chrA=10000), outfname)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

# And again...

bamFile <- regen(1000, c(chrA=10000), outfname)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

# Repeating with more reads.

bamFile <- regen(2000, c(chrA=10000), outfname)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

# Trying it out with multiple chromosomes.

bamFile <- regen(5000, c(chrA=10000, chrB=5000), outfname)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

# And again, with more reads.

bamFile <- regen(10000, c(chrA=10000, chrB=5000), outfname)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

# Trying it out with multiple BAM files.

bamFiles <- c(regen(500, c(chrA=1000, chrB=500), outfname), regen(500, c(chrA=1000, chrB=500), outfname))
comp(bamFiles, 50, cross=FALSE)
comp(bamFiles, 100, cross=FALSE)

# And again, with more reads.

bamFiles <- c(regen(5000, c(chrA=10000, chrB=5000), outfname), regen(5000, c(chrA=10000, chrB=5000), outfname))
comp(bamFiles, 50, cross=FALSE)
comp(bamFiles, 100, cross=FALSE)

###################################################################################################
# Throwing in some stress tests.

set.seed(789325)

# Where distance exceeds chromosome length.	

bamFile <- regen(10, c(chrA=20), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

bamFile <- regen(10, c(chrA=50), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

bamFile <- regen(5, c(chrA=100), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

# When the number of reads is zero in one chromosome.

bamFile <- regen(1, c(chrA=100, chrB=200), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

bamFile <- regen(2, c(chrA=100, chrB=200, chrC=30), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

bamFile <- regen(1, c(chrA=100, chrB=200, chrC=30), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

# Where the chromosome length is literally too small to compute correlations.

bamFile <- regen(1, c(chrA=1, chrB=1), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

bamFile <- regen(10, c(chrA=1, chrB=1), outfname)
comp(bamFile, 50, cross=TRUE)
comp(bamFile, 50, cross=FALSE)
comp(bamFile, 100, cross=FALSE)

###################################################################################################
# Cleaning out the directory.

unlink(fdir, recursive=TRUE);

###################################################################################################
# End.
