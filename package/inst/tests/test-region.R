###################################################################################################
# This script tests the region-counting capabilities of the 'csaw' package.

suppressWarnings(suppressPackageStartupMessages(library(csaw)))
source("simsam.R")

comp <- function(bamFiles, fraglen=200, right=0, left=0, spacing=20, filter=5, discard=GRanges(), restrict=NULL, forward=NA, final.len=NA) {
	if (length(final.len) && is.na(final.len)) { 
		ext <- rep(fraglen, length.out=length(bamFiles))
	} else {
		ext <- makeExtVector(fraglen, final.len)
	}

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
	
		# We compare windowCounts and regionCounts directly.
		repar <- readParam(discard=discard, restrict=restrict, minq=minq, dedup=dedup, forward=forward)
		x<-windowCounts(bamFiles, ext=ext, width=right+left+1, shift=left, spacing=spacing, 
			filter=filter, param=repar)
		all.regs <- rowRanges(x)
		if (!is.na(forward)) { strand(all.regs) <- "*" }

		y <- regionCounts(bamFiles, regions=all.regs, ext=ext, param=repar)
		if (!identical(assay(y), assay(x))) { stop("mismatch in count matrices") }
		if (!identical(y$totals, x$totals)) { stop("mismatch in total counts") }

		# If there's no rescaling, we pick a region in the middle and we check it with extractReads.
		chosen <- round(nrow(x)/2)
		my.reg <- all.regs[chosen]
	
		for (f in 1:length(bamFiles)) {
			collected <- extractReads(bamFiles[f], my.reg, param=repar, ext=ext[f])
			strand(collected) <- "*"
			if (!identical(assay(x)[chosen,f], length(collected))) { 
				stop("mismatch in the number of counts from extractReads")
			}
		}
	}

	return(head(assay(y)))
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

# More complex with variable fragment lengths.

comp(bamFiles, fraglen=c(100, 200), spacing=50)
comp(bamFiles, fraglen=c(100, 200), spacing=50, final.len=NULL)
comp(bamFiles, fraglen=c(100, 200), spacing=50, final.len=100)

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

comp(bamFiles, fraglen=c(100, 200), spacing=50)
comp(bamFiles, fraglen=c(100, 200), spacing=50, final.len=NULL)
comp(bamFiles, fraglen=c(100, 200), spacing=50, final.len=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, right=100)
comp(bamFiles, fraglen=200, left=10)
comp(bamFiles, fraglen=200, right=-5, left=10, spacing=20)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=20)
comp(bamFiles, fraglen=200, filter=40)

# One more time; sparse across the genome, but three files.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)
bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "C")))
comp(bamFiles, fraglen=100)
comp(bamFiles, fraglen=200)

comp(bamFiles, fraglen=c(100, 200, 150), spacing=50)
comp(bamFiles, fraglen=c(100, 200, 150), spacing=50, final.len=NULL)
comp(bamFiles, fraglen=c(100, 200, 150), spacing=50, final.len=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "C")))
comp(bamFiles, fraglen=100, left=50, spacing=100)
comp(bamFiles, fraglen=200, right=100, spacing=100)

bamFiles<-c(regen(3000, chromos, file.path(dir, "A")), regen(3000, chromos, file.path(dir, "B")), regen(3000, chromos, file.path(dir, "C")))
comp(bamFiles, fraglen=200, filter=10)
comp(bamFiles, fraglen=200, right=50, filter=50)

###################################################################################################
# Restricted and/or discarded.

chromos<-c(chrA=5000, chrB=5000, chrC=8000)

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, discard=makeDiscard(10, 200, chromos))
comp(bamFiles, fraglen=200, discard=makeDiscard(20, 100, chromos), restrict="chrA")
comp(bamFiles, fraglen=200, right=50, discard=makeDiscard(10, 200, chromos), restrict=c("chrA", "chrB"))

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=200, left=25, spacing=50, discard=makeDiscard(20, 200, chromos))
comp(bamFiles, fraglen=200, filter=1, discard=makeDiscard(5, 1000, chromos), restrict=c("chrC", "chrA"))
comp(bamFiles, fraglen=200, right=50, filter=2, discard=makeDiscard(20, 100, chromos))

bamFiles<-c(regen(100, chromos, file.path(dir, "A")), regen(100, chromos, file.path(dir, "B")))
comp(bamFiles, fraglen=100, filter=2, forward=TRUE)
comp(bamFiles, fraglen=100, filter=2, forward=FALSE, discard=makeDiscard(10, 200, chromos))
comp(bamFiles, fraglen=200, filter=2, left=20, spacing=50, forward=TRUE)
comp(bamFiles, fraglen=200, filter=2, left=20, spacing=50, forward=FALSE, discard=makeDiscard(20, 200, chromos))

###################################################################################################
# Cleaning up.

unlink(dir, recursive=TRUE)

###################################################################################################
# End.

