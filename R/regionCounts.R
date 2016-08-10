regionCounts <- function(bam.files, regions, ext=100, param=readParam())
# This just counts reads over regions. The only reason I'm using this and not
# some other package, is because (a) I want to avoid loading in more packages
# than I need, and (b) I need to count using the same reads (i.e., same values
# for 'ext', 'pe', and so on).
#
# written by Aaron Lun
# created 14 May 2014
# last modified 8 August 2016
{
	nbam <- length(bam.files)
    if (is.list(param)) { 
        .Deprecated(msg="supplying a list of readParam objects is deprecated, using first element only")
        param <- param[[1]]
    }
	extracted.chrs <- .activeChrs(bam.files, param$restrict)
	ext.data <- .collateExt(nbam, ext) 

	totals <- integer(nbam)
	nx <- length(regions)
	counts <- matrix(0L, nrow=nx, ncol=nbam)
	indices <- split(seq_len(nx), seqnames(regions))
    all.extras <- rep(list(list()), nbam)
    
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
        bp.out <- bpmapply(FUN=.region_counts, bam.file=bam.files, init.ext=ext.data$ext, 
                           MoreArgs=list(where=where, param=param, 
                                         final.ext=ext.data$final, outlen=outlen, 
                                         regions=regions, chosen=chosen),
                           BPPARAM=param$BPPARAM, SIMPLIFY=FALSE)

        for (bf in seq_along(bp.out)) {
            counts[chosen, bf] <- bp.out[[bf]]$counts
            totals[bf] <- totals[bf] + bp.out[[bf]]$totals
            all.extras[[bf]][[chr]] <- bp.out[[bf]]$extra
        }
	}

    strand(regions) <- .decideStrand(param)
	return(SummarizedExperiment(assays=SimpleList(counts=counts), 
		rowRanges=regions, 
		colData=.formatColData(bam.files, totals, ext.data, all.extras, param),
		metadata=list(final.ext=ext.data$final, param=param)))
}

.region_counts <- function(bam.file, where, param, 
                           init.ext, final.ext, outlen, 
                           regions, chosen) {
    if (param$pe!="both") {
        reads <- .getSingleEnd(bam.file, where=where, param=param)
        extended <- .extendSE(reads, ext=init.ext, final=final.ext, chrlen=outlen)
        frag.start <- extended$start
        frag.end <- extended$end

        extra <- cbind(c(mean(reads$forward$qwidth), mean(reads$reverse$qwidth)),
                       c(length(reads$forward$qwidth), length(reads$reverse$qwidth)))
    } else {
        out <- .getPairedEnd(bam.file, where=where, param=param)
        extra <- c(mean(out$size), length(out$size))

        checked <- .coerceFragments(out$pos, out$pos+out$size-1L, final=final.ext, chrlen=outlen)
        frag.start <- checked$start
        frag.end <- checked$end
    }
    
    # Counting the number of overlaps of any type with the known regions.
    if (length(chosen)==0L) { 
        counts <- 0L
    } else {
        counts <- countOverlaps(ranges(regions[chosen]), IRanges(frag.start, frag.end))
    }
    return(list(counts=counts, totals=length(frag.start), extra=extra))
}
