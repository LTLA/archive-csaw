windowCounts<-function(bam.files, spacing=50, left=0, right=0, ext=100, 
	pet=c("none", "both", "first", "second"), max.frag=500, filter=NULL, bin=NULL, 
	dedup=FALSE, minq=0, restrict=NULL, discard=NULL)
# Gets counts from BAM files at each position of the sliding window. Applies
# a gentle filter to remove the bulk of window positions with low counts.
# Returns a DGEList with count and total information, as well as a GRanges
# object specifying the window intervals.
# 
# written by Aaron Lun
# ages ago.
{   
	nbam<-length(bam.files)
	chrs<-scanBamHeader(bam.files[1])[[1]][[1]]
    if (!is.null(restrict)) { chrs<-chrs[names(chrs) %in% restrict] }
	discard <- .processDiscard(discard)

	# Processing input parameters.
	if (is.null(bin)) { 
		spacing<-as.integer(spacing+0.5)
		left<-as.integer(left+0.5)
		right<-as.integer(right+0.5)
		ext<-as.integer(ext+0.5)
		if (is.null(filter)) { filter <-5L*nbam }
	} else {
		# A convenience flag, which assigns sensible arguments to everything else.
		bin<-as.integer(bin+0.5)
		spacing<-bin
		right<-bin-1L
		left<-0L
		filter<-1L
		ext<-1L
	}
	pet <- match.arg(pet)
	max.frag <- as.integer(max.frag)
	minq <- as.integer(minq)
	dedup <- as.logical(dedup)

	# Initializing various collectable containers (non-empty so it'll work if no chromosomes are around).
	totals<-integer(nbam)	
	all.out<-list(matrix(0L, ncol=nbam, nrow=0))
	all.regions<-list(GRanges())
	ix<-1

	for (chr in names(chrs)) {
		outlen<-chrs[chr]		
		total.pts<-1L+as.integer((outlen-1L)/spacing)
		outcome<-matrix(0L, total.pts, nbam) 

		for (bf in 1:nbam) {
			where<-GRanges(chr, IRanges(1, outlen))
			if (pet!="both") {
				if (pet=="none") { 
   					reads<-.extractSET(bam.files[bf], where=where, dedup=dedup, minq=minq, discard=discard[[chr]])
				} else {
					reads <- .extractBrokenPET(bam.files[bf], where=where, dedup=dedup, minq=minq, discard=discard[[chr]], use.first=(pet=="first"))
				}
				frag.start<-ifelse(reads$strand=="+", reads$pos, reads$pos+reads$qwidth-ext)
				if (length(frag.start)) { frag.start<-pmin(frag.start, outlen) }
				frag.end<-frag.start+ext-1L
			} else {
				out<-.extractPET(bam.files[bf], where=where, dedup=dedup, minq=minq, discard=discard[[chr]])
				keep<-out$size <= max.frag
				frag.start<-out$pos[keep]

				# Only want to record each pair once in a bin, so forcing it to only use the starting read.
				if (!is.null(bin)) { frag.end<-frag.start }
				else { frag.end<-frag.start+out$size[keep]-1L; }
			}

			# Extending reads to account for window sizes > 1 bp. 'left' and 'right' refer to the length of
			# the window on either side of the center, so the start of each read must be extended by 'right'
			# and the end of each read must be extended by 'left'. We then pull out counts at the specified spacing.
			out<-.Call("R_get_rle_counts", frag.start-right, frag.end+left, total.pts, spacing, PACKAGE="csaw")
			if (is.character(out)) { stop(out) }
			outcome[,bf]<-out
			totals[bf]<-totals[bf]+length(frag.start)
		}

		# Filtering on row sums (for memory efficiency) and removing redundant windows. 
		center <- (1:total.pts-1L)*spacing+1L
		reg.start <- pmax(1L, center-left)
		reg.end <- pmin(outlen, center+right)
		keep <- rowSums(outcome)>=filter & c(TRUE, diff(reg.start)!=0L | diff(reg.end)!=0L)
		if (!sum(keep)) { next }
		all.out[[ix]]<-outcome[keep,,drop=FALSE]
		all.regions[[ix]]<-GRanges(chr, IRanges(reg.start[keep], reg.end[keep]))
		ix<-ix+1L
	}

	# Generating the remaining GRanges for output (suppressing numerous warnings).
	all.regions<-suppressWarnings(do.call(c, all.regions))
	seqlevels(all.regions)<-names(chrs)
	seqinfo(all.regions)<-Seqinfo(names(chrs), chrs)
	return(list(counts=do.call(rbind, all.out), totals=totals, region=all.regions))
}

########################################################

countWindows <- function(param, ...) 
# This is a wrapper for windowCounts which accepts a named parameter list and other
# non-default arguments. The latter will overwrite the former, if there are any 
# conflicts. The idea is to improve the re-callability of windowCounts.
#
# written by Aaron Lun
# 8 December, 2013
{
	other <- list(...)
	for (x in names(other)) { param[[x]]<-other[[x]] }
	do.call(windowCounts, param)
}

########################################################

.extractSET <- function(bam, where, dedup, minq, na.rm=TRUE, discard=NULL, ...) 
# Extracts single-end read data from a BAM file.
{
	reads<-scanBam(bam, param=ScanBamParam(what=c("strand", "pos", "qwidth", "mapq"),
			which=where, flag=scanBamFlag(isUnmappedQuery=FALSE, 
			isDuplicate=ifelse(dedup, FALSE, NA), ...)))[[1]]
    keep<-reads$mapq >= minq 
	if (na.rm) { keep<-keep & !is.na(reads$mapq) }
	reads$mapq <- NULL
	for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
	reads <- .discardReads(reads, discard)
	return(reads)
}

.processDiscard <- function(discard) 
# Processes a GRanges list for discarding objects.
{ 
	if (is.null(discard) || length(discard)==0L) { return(NULL) } 
	return(	split(ranges(discard), seqnames(discard)) )
}

.discardReads <- function(reads, discard.x) 
# Takes reads and discards them if they overlap with an intrachromosomal IRanges in 'discard.x'
{
	if (is.null(discard.x)) { return(reads) }
	keep <- !overlapsAny(IRanges(reads$pos, reads$pos+reads$qwidth-1L), discard.x)
	for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
	return(reads)
}

########################################################
