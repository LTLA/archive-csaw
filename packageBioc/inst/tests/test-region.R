###################################################################################################
# This script tests the region-counting capabilities of the 'csaw' package.

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

# We compare windowCounts and regionCounts directly.

comp <- function(bamFiles, fraglen=200, right=0, left=0, spacing=20, filter=5, discard=NULL, restrict=NULL) {
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
		y <- regionCounts(bamFiles, regions=x$region, ext=fraglen, discard=discard, restrict=restrict, minq=minq, dedup=dedup)
		if (!identical(y$counts, x$counts)) { stop("mismatch in count matrices") }
		if (!identical(y$totals, x$totals)) { stop("mismatch in total counts") }
	}
	return(head(y$counts))
}

###################################################################################################
# Setting up some variables to do the comparison.

dir<-"reg-test";
dir.create(dir);

set.seed(2123)
chromos<-c(chrA=10000, chrB=5000)

# Vanilla comparison.

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, spacing=20)
comp(bamFiles, fraglen=200, spacing=50)

# More complex with right arguments.

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, right=30, spacing=20)
comp(bamFiles, fraglen=200, left=5, spacing=25)
comp(bamFiles, fraglen=150, right=-10, left=10, spacing=30)

# Even more complex, with filtering arguments

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=10)
comp(bamFiles, fraglen=200, filter=15)
comp(bamFiles, fraglen=200, filter=20)

# And again, with a different chromosome set-up.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, left=10)
comp(bamFiles, fraglen=200, right=-5, left=10, spacing=20)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=20)
comp(bamFiles, fraglen=200, filter=40)

# One more time; sparse across the genome, but three files.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, left=50, spacing=100)
comp(bamFiles, fraglen=200, right=100, spacing=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=200, filter=10)
comp(bamFiles, fraglen=200, right=50, filter=50)

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
comp(bamFiles, fraglen=200, right=50, discard=makeDiscard(10, 200), restrict=c("chrA", "chrB"))

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=200, left=25, spacing=50, discard=makeDiscard(20, 200))
comp(bamFiles, fraglen=200, filter=1, discard=makeDiscard(5, 1000), restrict=c("chrC", "chrA"))
comp(bamFiles, fraglen=200, right=50, filter=2, discard=makeDiscard(20, 100))

###################################################################################################
# Cleaning up.

unlink(dir, recursive=TRUE)

###################################################################################################
# End.

