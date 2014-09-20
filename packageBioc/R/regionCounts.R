regionCounts <- function(bam.files, regions, ext=100, param=readParam())
# This just counts reads over regions. The only reason I'm using this and not
# some other package, is because (a) I want to avoid loading in more packages
# than I need, and (b) I need to count using the same reads (i.e., same values
# for ext, pet, and so on).
#
# written by Aaron Lun
# 14th May 2014
{
    nbam <- length(bam.files)
	extracted <- .processIncoming(bam.files, param$restrict, param$discard)
    pet <- param$pet
    max.frag <- param$max.frag
	minq <- param$minq
	dedup <- param$dedup
	rescue.pairs <- param$rescue.pairs
	rescue.ext <- param$rescue.ext
	ext <- as.integer(ext)

    totals <- integer(nbam)
	nx <- length(regions)
	counts <- matrix(0L, nrow=nx, ncol=nbam)
	indices <- split(1:nx, seqnames(regions))

    for (chr in names(extracted$chrs)) {
		chosen <- indices[[chr]]
		if (length(chosen)==0L) { next }
        outlen <- extracted$chrs[[chr]]
        where <- GRanges(chr, IRanges(1, outlen))

		# Pulling out reads as previously described.
        for (bf in 1:nbam) {
            if (pet!="both") {
                if (pet=="none") {
                    reads <- .extractSET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        discard=extracted$discard[[chr]])
                } else {
                    reads <- .extractBrokenPET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        discard=extracted$discard[[chr]], use.first=(pet=="first"))
                }
                frag.start <- ifelse(reads$strand=="+", reads$pos, reads$pos+reads$qwidth-ext)
                if (length(frag.start)) { frag.start <- pmin(frag.start, outlen) }
                frag.end <- frag.start+ext-1L
            } else {
                if (rescue.pairs) {
                    out <- .rescuePET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        max.frag=max.frag, ext=rescue.ext, discard=extracted$discard[[chr]])
                } else {
                    out <- .extractPET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        discard=extracted$discard[[chr]], max.frag=max.frag)
                }
                frag.start <- out$pos
				frag.end <- frag.start+out$size-1L
            }
		
			# Counting the number of overlaps of any type with the known regions.
			totals[bf] <- totals[bf] + length(frag.start)
			counts[chosen,bf] <- countOverlaps(ranges(regions[chosen]), IRanges(frag.start, frag.end))
		}
	}
	return(SummarizedExperiment(assays=counts, 
			rowData=regions, colData=DataFrame(totals=totals)))
}
