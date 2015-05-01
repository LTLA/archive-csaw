checkBimodality <- function(bam.files, regions, width=100, param=readParam(), prior.count=2) 
# This gives the maximum strand bimodality score for any base pair in each
# region. The idea is to try to distinguish between genuine TF binding sites
# and read stacks or other artifacts.
#
# written by Aaron Lun
# created 1 May 2015
{
	nbam <- length(bam.files)
	paramlist <- .makeParamList(nbam, param)
	extracted.chrs <- .activeChrs(bam.files, paramlist[[1]]$restrict)
	ext.data <- .collateExt(nbam, width) 

	totals <- integer(nbam)
	nx <- length(regions)
	out.scores <- rep(NA, nx)
	indices <- split(1:nx, seqnames(regions))
	prior.count <- as.double(prior.count)

	for (chr in names(extracted.chrs)) {
		chosen <- indices[[chr]]
		outlen <- extracted.chrs[[chr]]
		where <- GRanges(chr, IRanges(1, outlen))

		# Pulling out reads as previously described.
		collected <- list()
		for (bf in 1:nbam) {
			curpar <- paramlist[[bf]]
			if (curpar$pe=="both") { 
				stop("bimodality checking not supported for paired-end mode") 
			} else if (curpar$pe=="none") {
				reads <- .extractSE(bam.files[bf], where=where, param=curpar)
			} else {
				reads <- .extractBrokenPE(bam.files[bf], where=where, param=curpar)
			}
			
			is.reverse <- reads$strand!="+"
			five.prime <- reads$pos
			five.prime[is.reverse] <- reads$pos[is.reverse] + reads$qwidth[is.reverse] - 1L

			curwidth <- ext.data$ext[bf]
			start.point <- five.prime - curwidth + 1L
			if (curwidth < 1L || !is.finite(curwidth)) { stop('width must be a non-negative integer') }

			o <- order(start.point)
			collected[[bf]] <- list(start.point[o], as.integer(!is.reverse[o]), curwidth)
		}

		# Computing bimodality scores.
		rstarts <- start(regions)[chosen]
		rends <- end(regions)[chosen]
		ro <- order(rstarts)
		out <- .Call(cxx_check_bimodality, collected, rstarts[ro], rends[ro], prior.count)
		if (is.character(out)) { stop(out) }

		out.scores[chosen][ro] <- out 
	}

	return(out.scores)	
}
