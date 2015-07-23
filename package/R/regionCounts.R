regionCounts <- function(bam.files, regions, ext=100, param=readParam())
# This just counts reads over regions. The only reason I'm using this and not
# some other package, is because (a) I want to avoid loading in more packages
# than I need, and (b) I need to count using the same reads (i.e., same values
# for 'ext', 'pe', and so on).
#
# written by Aaron Lun
# created 14 May 2014
# last modified 22 July 2015
{
	nbam <- length(bam.files)
	paramlist <- .makeParamList(nbam, param)
	extracted.chrs <- .activeChrs(bam.files, paramlist[[1]]$restrict)
	ext.data <- .collateExt(nbam, ext) 

	totals <- integer(nbam)
	nx <- length(regions)
	counts <- matrix(0L, nrow=nx, ncol=nbam)
	indices <- split(seq_len(nx), seqnames(regions))

	# No sense in doing so; you can set param$forward for strand-specific counting.
	if (any(strand(regions)!="*")) { 
		warning("ignoring strandedness of supplied regions") 
		strand(regions) <- "*"
	}

	for (chr in names(extracted.chrs)) {
		chosen <- indices[[chr]]
		outlen <- extracted.chrs[[chr]]
		where <- GRanges(chr, IRanges(1, outlen))

		# Pulling out reads as previously described.
		for (bf in seq_len(nbam)) {
			curpar <- paramlist[[bf]]
			if (curpar$pe!="both") {
				reads <- .getSingleEnd(bam.files[bf], where=where, param=curpar)
				extended <- .extendSE(reads, ext=ext.data$ext[bf], final=ext.data$final, chrlen=outlen)
				frag.start <- extended$start
				frag.end <- extended$end
			} else {
				out <- .getPairedEnd(bam.files[bf], where=where, param=curpar)
				checked <- .checkFragments(out$pos, out$pos+out$size-1L, final=ext.data$final, chrlen=outlen)
   				frag.start <- checked$start
				frag.end <- checked$end
			}
		
			# Counting the number of overlaps of any type with the known regions.
			totals[bf] <- totals[bf] + length(frag.start)
			if (length(chosen)==0L) { next }
			counts[chosen,bf] <- countOverlaps(ranges(regions[chosen]), IRanges(frag.start, frag.end))
		}
	}

	strand(regions) <- .decideStrand(paramlist)
	dim(paramlist) <- c(nbam, 1)
	colnames(paramlist) <- "param"
	return(SummarizedExperiment(assays=SimpleList(counts=counts), 
		rowRanges=regions, 
		colData=DataFrame(bam.files, totals=totals, ext=ext.data$ext, paramlist),
		metadata=list(final.ext=ext.data$final)))
}
