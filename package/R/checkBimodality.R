checkBimodality <- function(bam.files, regions, width=100, param=readParam(), 
	prior.count=2, invert=FALSE) 
# This gives the maximum strand bimodality score for any base pair in each
# region. The idea is to try to distinguish between genuine TF binding sites
# and read stacks or other artifacts.
#
# written by Aaron Lun
# created 1 May 2015
# last modified 22 July 2015
{
	nbam <- length(bam.files)
	paramlist <- .makeParamList(nbam, param)
	extracted.chrs <- .activeChrs(bam.files, paramlist[[1]]$restrict)
	ext.data <- .collateExt(nbam, width) 
	invert <- as.logical(invert)

	totals <- integer(nbam)
	nx <- length(regions)
	out.scores <- rep(NA_real_, nx)
	indices <- split(seq_len(nx), seqnames(regions))
	prior.count <- as.double(prior.count)

	for (chr in names(extracted.chrs)) {
		chosen <- indices[[chr]]
		if (length(chosen)==0L) { next } 
		outlen <- extracted.chrs[[chr]]
		where <- GRanges(chr, IRanges(1, outlen))

		# Pulling out reads as previously described.
		collected <- list()
		for (bf in seq_len(nbam)) {
			curpar <- paramlist[[bf]]
    
       		if (curpar$pe=="both") {
				out <- .getPairedEnd(bam.files[bf], where=where, param=curpar, with.reads=TRUE)
				if (.needsRescue(curpar)) { 
					reads <- mapply(c, out$left, out$right, out$rescued, SIMPLIFY=FALSE)
				} else {
					reads <- mapply(c, out$left, out$right, SIMPLIFY=FALSE)
				}
			} else {
				reads <- .getSingleEnd(bam.files[bf], where=where, param=curpar)
			}
			
			is.forward <- as.integer(reads$strand=="+")
			left.pos <- reads$pos
			right.pos <- reads$pos + reads$qwidth - 1L 

			if (is.unsorted(left.pos)) { 
				o <- order(left.pos)
				left.pos <- left.pos[o]
				right.pos <- right.pos[o]
				is.forward <- is.forward[o]
			}

			curwidth <- ext.data$ext[bf]
			if (curwidth < 1L || !is.finite(curwidth)) { stop('width must be a non-negative integer') }
			collected[[bf]] <- list(left.pos, right.pos, is.forward, curwidth)
		}

		# Checking region order.
		rstarts <- start(regions)[chosen]
		rends <- end(regions)[chosen]
		ro <- order(rstarts)
		
		# Computing bimodality scores.
		out <- .Call(cxx_check_bimodality, collected, rstarts[ro], rends[ro], prior.count, invert)
		if (is.character(out)) { stop(out) }
		out.scores[chosen][ro] <- out 
	}

	return(out.scores)	
}
