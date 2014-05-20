windowCounts<-function(bam.files, spacing=50, width=1, ext=100, shift=0,
	pet=c("none", "both", "first", "second"), max.frag=500, rescue.pairs=FALSE,
	filter=NULL, bin=FALSE, dedup=FALSE, minq=0, restrict=NULL, discard=NULL)
# Gets counts from BAM files at each position of the sliding window. Applies
# a gentle filter to remove the bulk of window positions with low counts.
# Returns a DGEList with count and total information, as well as a GRanges
# object specifying the window intervals.
# 
# written by Aaron Lun
# ages ago.
{   
	nbam<-length(bam.files)
	extracted <- .processIncoming(bam.files, restrict, discard)

	# Processing input parameters.
	if (length(bin)>1L || !is.logical(bin)) { stop("bin must be a logical scalar") }
	if (width < 1) { stop("width must be a positive integer") }
	if (!bin) { 
		spacing <- as.integer(spacing+0.5)
		left <- as.integer(shift+0.5)
		right <- as.integer(width+0.5) - left - 1L
		ext <- as.integer(ext+0.5)
		if (is.null(filter)) { filter <- 5L*nbam }
	} else {
		# A convenience flag, which assigns sensible arguments to everything else.
		spacing <- as.integer(width+0.5)
		left <- as.integer(shift+0.5)
		right <- spacing - 1L - left
		filter <- ext <- 1L
	}
	pet <- match.arg(pet)
	max.frag <- as.integer(max.frag)
	minq <- as.integer(minq)
	dedup <- as.logical(dedup)

# Figuring out what the extensions are necessary. We've reparameterised it so
# that 'left' and 'right' refer to the extension of the window from a nominal
# 'center' point. This simplifies counting. However, we do need to account for
# loss of a point from the front when shift/left is non-zero. There's also a 
# gain but that's dealt with later.
	if (left >= spacing) { stop("shift must be less than the spacing") }
	if (left < 0L) { stop("shift must be positive") }
	at.start <- right >= 0L
	first.pt <- ifelse(at.start, 1L, spacing+1L)

	# Initializing various collectable containers (non-empty so it'll work if no chromosomes are around).
	totals<-integer(nbam)	
	all.out<-list(matrix(0L, ncol=nbam, nrow=0))
	all.regions<-list(GRanges())
	ix<-1

	for (chr in names(extracted$chrs)) {
		outlen<-extracted$chrs[[chr]]		
		where<-GRanges(chr, IRanges(1, outlen))
		# Accounting for the gain of a point from the back when shift/left is non-zero. 
		at.end <- spacing - (outlen - 1L) %% spacing <= left
		total.pts <- as.integer((outlen-1)/spacing) + at.start + at.end
		outcome <- matrix(0L, total.pts, nbam) 

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

				# Only want to record each pair once in a bin, so forcing it to only use the 5' end.
				if (bin) { frag.end<-frag.start }
				else { frag.end<-frag.start+out$size-1L }
			}

# Extending reads to account for window sizes > 1 bp. The start of each read
# must be extended by 'right' and the end of each read must be extended by
# 'left'. We then pull out counts at the specified spacing. We do have to 
# keep track of whether or not we want to use the first point, though.
			out<-.Call("R_get_rle_counts", frag.start-right, frag.end+left, total.pts, spacing, at.start, PACKAGE="csaw")
			if (is.character(out)) { stop(out) }
			outcome[,bf]<-out
			totals[bf]<-totals[bf]+length(frag.start)
		}

# Filtering on row sums (for memory efficiency). Note that redundant windows
# are avoided by enforcing 'shift < spacing', as a window that is shifted to 
# cover the whole chromosome cannot be spaced to a new position where it still 
# covers the whole chromosome. The next one must be inside the chromosome.
		keep <- rowSums(outcome)>=filter 
		if (!any(keep)) { next } 
		else if (!all(keep)) { outcome <- outcome[keep,,drop=FALSE] }
		all.out[[ix]] <- outcome

		center <- (which(keep)-1L)*spacing + first.pt
		reg.start <- pmax(1L, center-left)
		reg.end <- pmin(outlen, center+right)
		all.regions[[ix]]<-GRanges(chr, IRanges(reg.start, reg.end))
		ix<-ix+1L
	}

	# Generating the remaining GRanges for output (suppressing numerous warnings).
	all.regions<-suppressWarnings(do.call(c, all.regions))
	seqlevels(all.regions) <- names(extracted$chrs)
	seqlengths(all.regions) <- extracted$chrs
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

.extractSET <- function(bam, where, dedup, minq, discard=NULL, extras="strand", ...) 
# Extracts single-end read data from a BAM file with removal of unmapped,
# duplicate and poorly mapped/non-unique reads. We also discard reads in the
# specified discard regions. In such cases, the offending reads must be wholly
# within the repeat region.  We use the real alignment width, just in case we
# have very long reads in the alignment that are heavily soft-clipped (i.e., they
# should be reported as within but the read length will put them out).
{
	all.fields <- c("pos", "qwidth", "mapq", extras)
	if (!is.null(discard)) { all.fields <- c(all.fields, "cigar") }	
	all.fields <- unique(all.fields)
	reads <- scanBam(bam, param=ScanBamParam(what=all.fields,
		which=where, flag=scanBamFlag(isUnmappedQuery=FALSE, 
		isDuplicate=ifelse(dedup, FALSE, NA), ...)))[[1]]
   
	# Filtering by MAPQ.
	keep<-reads$mapq >= minq & !is.na(reads$mapq) 
	if (!"mapq" %in% extras) { reads$mapq <- NULL }
	for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
	
	# Filtering by discard regions. Using alignment width so long reads can escape repeats.
	if (!is.null(discard)) {
 	   	require(GenomicAlignments)	
		awidth <- cigarWidthAlongReferenceSpace(reads$cigar)
		keep <- !overlapsAny(IRanges(reads$pos, reads$pos+awidth-1L), discard, type="within")
		for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
		if (!"cigar" %in% extras) { reads$cigar <- NULL }
	}
	return(reads)
}

.processIncoming <- function(bam.files, restrict, discard) 
# Processes the incoming data; checks that bam headers are all correct,
# truncates the list according to 'restrict', processes the GRanges list 
# for discarding objects.
{ 
	originals <- NULL
	for (bam in bam.files) {
		chrs <- scanBamHeader(bam)[[1]][[1]]
		chrs <- chrs[order(names(chrs))]
		if (is.null(originals)) { originals <- chrs } 
		else if (!identical(originals, chrs)) { 
			warning("chromosomes are not identical between BAM files")
			pairing <- match(names(originals), names(chrs))
			originals <- pmin(originals[!is.na(pairing)], chrs[pairing[!is.na(pairing)]])
		}
	}

	if (!is.null(restrict)) { originals <- originals[names(originals) %in% restrict] }
	if (!is.null(discard)) {
		discard <- split(ranges(discard), seqnames(discard), drop=TRUE)
		for (x in names(discard)) { 
			if (x %in% names(originals)) { 
				discard[[x]] <- reduce(discard[[x]])
			} else {
				discard[[x]] <- NULL 
			}
		}
	} 
	return(list(discard=discard, chrs=originals))
}

########################################################
