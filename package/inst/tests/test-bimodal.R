# This checks the checkBimodality function in the csaw package.

suppressWarnings(suppressPackageStartupMessages(library(csaw)))
source("simsam.R")

comp <- function(bamFiles, regions, width=100, discard=GRanges(), restrict=NULL, prior.count=2) { 
	width <- rep(width, length.out=length(bamFiles))
	by.chr <- split(1:length(regions), seqnames(regions))
	chr.sizes <- Rsamtools::scanBamHeader(bamFiles[1])[[1]][[1]]
	chr.in.use <- names(by.chr)
	if (!is.null(restrict)) { chr.in.use <- intersect(chr.in.use, restrict) } 

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
        repar <- readParam(discard=discard, restrict=restrict, minq=minq, dedup=dedup)
		observed <- checkBimodality(bamFiles, regions, param=repar, width=width, prior.count=prior.count)
		output <- rep(NA, length(regions))

		for (chr in chr.in.use) { 
			full.chr <- GRanges(chr, IRanges(1, chr.sizes[[chr]]))
			left.forward <- right.forward <- left.reverse <- right.reverse <- 0

			for (b in 1:length(bamFiles)) { 
				all.reads <- extractReads(bamFiles[b], full.chr, param=repar)

				# Computing coverage, somewhat counterintuitively; 'left.forward', for example, is extended to 
				# the right, as it must find all overlaps to the left of the region of interest.
				is.forward <- strand(all.reads)=="+"
				forward.reads <- all.reads[is.forward]
				reverse.reads <- all.reads[!is.forward]
	
				start.F <- start(forward.reads)
				end.F <- end(forward.reads)
				left.forward <- left.forward + coverage(IRanges(start.F, end.F+width[b]-1), width=end(full.chr))
				right.forward <- right.forward + coverage(IRanges(start.F-width[b]+1, end.F), width=end(full.chr))

				end.R <- end(reverse.reads)
				start.R <- start(reverse.reads)
				left.reverse <- left.reverse + coverage(IRanges(start.R, end.R+width[b]-1), width=end(full.chr))
				right.reverse <- right.reverse + coverage(IRanges(start.R-width[b]+1, end.R), width=end(full.chr))
			}

			for (r in by.chr[[chr]]) { 
				cur.reg <- regions[r]

				# Computing the bimodality statistic.
				lf.set <- as.integer(left.forward[start(cur.reg):end(cur.reg)])
				rf.set <- as.integer(right.forward[start(cur.reg):end(cur.reg)])
				lr.set <- as.integer(left.reverse[start(cur.reg):end(cur.reg)])
				rr.set <- as.integer(right.reverse[start(cur.reg):end(cur.reg)])

	  			bistat <- max(pmin((lf.set+prior.count)/(lr.set+prior.count), (rr.set+prior.count)/(rf.set+prior.count)))
  				output[r] <- bistat
			}
		}

		# Comparing values.
		is.lost <- is.na(output)
		if (!identical(is.lost, is.na(observed))) { stop("mismatch in missing values") }
		if (any(abs(output[!is.lost] - observed[!is.lost]) > (0.001 + output[!is.lost])*1e-6)) { 
			stop("mismatch in computed bimodality values") 
		}
	}

	return(head(observed))
}

####################################################################################################

set.seed(3985)

dir <- "bimodal-test"
dir.create(dir)
chromos<-c(chrA=5000, chrB=5000, chrC=8000)

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
my.ranges <- generateWindows(chromos, 10, 10)
comp(bamFiles, my.ranges)
comp(bamFiles, my.ranges, width=200)
comp(bamFiles, my.ranges, width=c(100, 200))
comp(bamFiles, my.ranges, prior.count=5)
comp(bamFiles, my.ranges, discard=makeDiscard(10, 200, chromos))
comp(bamFiles, my.ranges, restrict=c("chrC", "chrB"))

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
my.ranges <- generateWindows(chromos, 10, 500)
comp(bamFiles, my.ranges)
comp(bamFiles, my.ranges, width=200)
comp(bamFiles, my.ranges, width=c(100, 200))
comp(bamFiles, my.ranges, prior.count=5)
comp(bamFiles, my.ranges, discard=makeDiscard(10, 200, chromos))
comp(bamFiles, my.ranges, restrict=c("chrC", "chrB"))

bamFiles<-c(regen(1000, chromos, file.path(dir, "A")), regen(1000, chromos, file.path(dir, "B")))
my.ranges <- generateWindows(chromos, 20, round(runif(20, 100, 1000)))
comp(bamFiles, my.ranges)
comp(bamFiles, my.ranges, width=200)
comp(bamFiles, my.ranges, width=c(100, 200))
comp(bamFiles, my.ranges, prior.count=5)
comp(bamFiles, my.ranges, discard=makeDiscard(10, 200, chromos))
comp(bamFiles, my.ranges, restrict=c("chrC", "chrB"))

####################################################################################################
# Cleaning up.

unlink(dir, recursive=TRUE)

####################################################################################################
# End.
