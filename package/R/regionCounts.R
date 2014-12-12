regionCounts <- function(bam.files, regions, ext=NULL, param=readParam())
# This just counts reads over regions. The only reason I'm using this and not
# some other package, is because (a) I want to avoid loading in more packages
# than I need, and (b) I need to count using the same reads (i.e., same values
# for ext, pet, and so on).
#
# written by Aaron Lun
# created 14 May 2014
# last modified 12 December 2014
{
	nbam <- length(bam.files)
	paramlist <- .makeParamList(nbam, param)
	if (!is.null(ext)) { paramlist <- reformList(paramlist, ext=ext) }
	extracted.chrs <- .activeChrs(bam.files, paramlist[[1]]$restrict)

    totals <- integer(nbam)
	nx <- length(regions)
	counts <- matrix(0L, nrow=nx, ncol=nbam)
	indices <- split(1:nx, seqnames(regions))

    for (chr in names(extracted.chrs)) {
		chosen <- indices[[chr]]
        outlen <- extracted.chrs[[chr]]
        where <- GRanges(chr, IRanges(1, outlen))

		# Pulling out reads as previously described.
        for (bf in 1:nbam) {
			curpar <- paramlist[[bf]]
            if (curpar$pet!="both") {
                if (curpar$pet=="none") {
                    reads <- .extractSET(bam.files[bf], where=where, param=curpar)
                } else {
                    reads <- .extractBrokenPET(bam.files[bf], where=where, param=curpar)
                }
				extended <- .extendSE(reads, chrlen=outlen, param=curpar)
				frag.start <- extended$start
				frag.end <- extended$end
            } else {
                if (curpar$rescue.pairs) {
                    out <- .rescuePET(bam.files[bf], where=where, param=curpar)
                } else {
                    out <- .extractPET(bam.files[bf], where=where, param=curpar)
                }
                frag.start <- out$pos
				frag.end <- frag.start+out$size-1L
            }
		
			# Counting the number of overlaps of any type with the known regions.
			totals[bf] <- totals[bf] + length(frag.start)
			if (length(chosen)==0L) { next }
			counts[chosen,bf] <- countOverlaps(ranges(regions[chosen]), IRanges(frag.start, frag.end))
		}
	}

	if (is.list(param)) { 
		index <- 1:nbam
	} else {
		index <- 1L
		paramlist <- paramlist[1]
	}
	return(SummarizedExperiment(assays=counts, 
		rowData=regions, 
		colData=DataFrame(bam.files, totals=totals, param=index),
		exptData=SimpleList(param=paramlist)))
}
