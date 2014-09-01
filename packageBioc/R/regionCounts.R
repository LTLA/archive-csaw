regionCounts <- function(bam.files, regions,  ext=100, 
    pet=c("none", "both", "first", "second"), max.frag=500, rescue.pairs=FALSE,
	dedup=FALSE, minq=NA, restrict=NULL, discard=NULL)
# This just counts reads over regions. The only reason I'm using this and not
# some other package, is because (a) I want to avoid loading in more packages
# than I need, and (b) I need to count using the same reads (i.e., same values
# for ext, pet, and so on).
#
# written by Aaron Lun
# 14th May 2014
{
    nbam<-length(bam.files)
	extracted <- .processIncoming(bam.files, restrict, discard)
    pet <- match.arg(pet)
    max.frag <- as.integer(max.frag)
	minq <- as.integer(minq)
	dedup <- as.logical(dedup)
	ext <- as.integer(ext)

    totals <- integer(nbam)
	nx <- length(regions)
	counts <- matrix(0L, nrow=nx, ncol=nbam)
	indices <- split(1:nx, seqnames(regions))

    for (chr in names(extracted$chrs)) {
		chosen <- indices[[chr]]
		if (is.null(chosen)) { next }
        outlen<-extracted$chrs[[chr]]
        where<-GRanges(chr, IRanges(1, outlen))

		# Pulling out reads as previously described.
        for (bf in 1:nbam) {
            if (pet!="both") {
                if (pet=="none") {
                    reads<-.extractSET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        discard=extracted$discard[[chr]])
                } else {
                    reads <- .extractBrokenPET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        discard=extracted$discard[[chr]], use.first=(pet=="first"))
                }
                frag.start<-ifelse(reads$strand=="+", reads$pos, reads$pos+reads$qwidth-ext)
                if (length(frag.start)) { frag.start<-pmin(frag.start, outlen) }
                frag.end<-frag.start+ext-1L
            } else {
                if (rescue.pairs) {
                    out<-.rescuePET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        max.frag=max.frag, ext=ext, discard=extracted$discard[[chr]])
                } else {
                    out<-.extractPET(bam.files[bf], where=where, dedup=dedup, minq=minq,
                        discard=extracted$discard[[chr]], max.frag=max.frag)
                }
                frag.start<-out$pos
				frag.end<-frag.start+out$size-1L
            }
		
			# Counting the number of overlaps of any type with the known regions.
			totals[bf] <- totals[bf] + length(frag.start)
			counts[chosen,bf] <- countOverlaps(ranges(regions[chosen]), IRanges(frag.start, frag.end))
		}
	}
	return(list(counts=counts, totals=totals))
}
