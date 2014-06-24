plotRegion <- function(cur.region, bam.file, dedup=FALSE, minq=NA, discard=NULL, 
	pet=c("none", "both", "first", "second"), max.frag=500,	rescue.ext=NULL,
	max.depth=NULL, fcol="blue", rcol="red", 
	xlab="Genomic position (bp)", ylab="Read depth", ...)
# Exactly as specified. Takes a region and plots it in bimodal mode, with
# options for duplicate removal, mapping quality enhancement, colorization,
# and PET manipulation.
#
# written by Aaron Lun
{
    if (length(cur.region)!=1L) { stop("exactly one range is required for plotting") }
    chrs<-scanBamHeader(bam.file)[[1]][[1]]
	cur.chr<-as.character(seqnames(cur.region)[1])
	if (!(cur.chr %in% names(chrs))) { stop("cannot find current chromosome in the BAM file header") }

	# Expanding the extracted region so that it's continuous past the plot margins.
	cur.width<-width(cur.region)
	max.len<-chrs[[cur.chr]]
	ext<-cur.width/2+1000 + ifelse(pet=="both", max.frag, 0) 
	actual.region<-GRanges(cur.chr, IRanges(pmax(1L, start(cur.region)-ext), 
		pmin(max.len, end(cur.region)+ext)))

	# Dropping additional reads if required.
    if (!is.null(discard)) { discard <- ranges(discard[overlapsAny(discard, actual.region)]) }

	# Pulling out reads from a region and setting up coverage RLE's.
	pet <- match.arg(pet)
	if (pet!="both") {
		if (pet=="none") { 
			cur.reads <- .extractSET(bam.file, where=actual.region, dedup=dedup, 
				minq=minq, discard=discard)
		} else {
			cur.reads <- .extractBrokenPET(bam.file, where=actual.region, dedup=dedup, 
				minq=minq, discard=discard, use.first=(pet=="first"))
		}
		forward<-cur.reads$strand=="+"
  		starts<-cur.reads$pos
		ends<-starts+cur.reads$qwidth-1L
	} else {
		if (!is.null(rescue.ext)) {
			cur.reads <- .rescuePET(bam.file, where=actual.region, dedup=dedup, 
					minq=minq, discard=discard, max.frag=max.frag, ext=rescue.ext)
		} else {
			cur.reads <- .extractPET(bam.file, where=actual.region, dedup=dedup, 
					minq=minq, discard=discard, max.frag=max.frag)
		}

		# Plotting PETs with a bit more care.
		nx <- length(cur.reads$pos)
		forward <- c(logical(nx), !logical(nx))
		midpoint <- as.integer(cur.reads$pos + cur.reads$size/2)
		starts <- c(midpoint+1L, cur.reads$pos)
 	    ends <- c(cur.reads$pos + cur.reads$size, midpoint)
	}
	pos<-coverage(IRanges(starts[forward], ends[forward]), width=max.len)
	neg<-coverage(IRanges(starts[!forward], ends[!forward]), width=max.len)
	
	# Just plotting.
	if (is.null(max.depth)) { max.depth<-max(runValue(pos), runValue(neg)) }
	ylim<-c(-max.depth, max.depth)
	plot(-1, 1, xlim=c(start(cur.region), end(cur.region)), ylim=ylim, type="n",
			xlab=xlab, ylab=ylab, ...)
	rect(start(pos), 0, c(start(pos)[-1], length(pos)), runValue(pos), col=fcol, border=NA)
	rect(start(neg), 0, c(start(neg)[-1], length(neg)), -runValue(neg), col=rcol, border=NA)

	return(invisible(NULL))
}

