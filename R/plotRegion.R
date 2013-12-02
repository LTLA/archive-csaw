plotRegion <- function(cur.region, bam.file, dedup=FALSE, minq=0, max.depth=NULL, fcol="blue", rcol="red", 
	xlab="Genomic position (bp)", ylab="Read depth", ...)
# Exactly as specified. Takes a region and plots it in bimodal mode, with options for duplicate 
# removal and colorization.
#
# written by Aaron Lun
{
    if (length(cur.region)!=1L) { stop("exactly one range is required for plotting") }
	cur.width<-width(cur.region)
	if (cur.width < 50) { 
		warning("width of specified region may be too low")
		if (cur.width==1) { 
			cur.width<-50
			end(cur.region)<-start(cur.region)+cur.width
		}
	}

    chrs<-scanBamHeader(bam.file)[[1]][[1]]
	cur.chr<-as.character(seqnames(cur.region)[1])
	if (!(cur.chr %in% names(chrs))) { stop("cannot find current chromosome in the BAM file header") }
	max.len<-chrs[[cur.chr]]
	ext<-cur.width/2+1000
	actual.region<-GRanges(cur.chr, IRanges(pmax(1L, start(cur.region)-ext), 
		pmin(max.len, end(cur.region)+ext)))

	# Pulling out reads from a region and setting up coverage RLE's.
	cur.reads<-.extractSET(bam.file, where=actual.region, dedup=dedup, minq=minq)
	forward<-cur.reads$strand=="+"
  	starts<-cur.reads$pos
	ends<-starts+cur.reads$qwidth-1L
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

