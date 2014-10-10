extractReads <- function(cur.region, bam.file, param=readParam())
# Exactly as specified. Takes a region and plots it in bimodal mode, with
# options for duplicate removal, mapping quality enhancement, colorization,
# and PET manipulation.
#
# written by Aaron Lun
# 1 September 2014
{
    if (length(cur.region)!=1L) { stop("exactly one range is required for plotting") }
    chrs <- scanBamHeader(bam.file)[[1]][[1]]
	cur.chr <- as.character(seqnames(cur.region)[1])
	if (length(param$restrict) && ! cur.chr %in% param$restrict) { stop("current chromosome not in restricted subset") }
	if (! cur.chr %in% names(chrs)) { stop("cannot find current chromosome in the BAM file header") }
	max.len <- chrs[[cur.chr]]
	sqi <- Seqinfo(cur.chr, max.len)

	# Extracting all-of-chromosome for paired-end rescue, as you need to find the read with the higher MAPQ.
	expand <- 0L
	if (param$pet=="both" && param$rescue.pairs) {
		actual.region <- GRanges(cur.chr, IRanges(1L, max.len)) 
	} else {
		if (param$pet=="both") {
			actual.region <- GRanges(cur.chr, IRanges(max(1L, start(cur.region)-param$max.frag),
				min(max.len, end(cur.region)+param$max.frag)))
		} else {
			actual.region <- cur.region
			seqlevels(actual.region) <- cur.chr
		}
	}

	# Dropping additional reads if required.
	minq <- param$minq
	dedup <- param$dedup
    if (length(param$discard)) { 
		cur.lost <- param$discard[seqnames(param$discard)==cur.chr]
		seqlevels(cur.lost) <- cur.chr
		discard <- ranges(cur.lost[overlapsAny(cur.lost, actual.region)]) 
	} else { discard <- NULL }

	# Pulling out reads from a region and setting up coverage RLE's.
	if (param$pet!="both") {
		if (param$pet=="none") { 
			cur.reads <- .extractSET(bam.file, where=actual.region, dedup=dedup, 
				minq=minq, discard=discard)
		} else {
			cur.reads <- .extractBrokenPET(bam.file, where=actual.region, dedup=dedup, 
				minq=minq, discard=discard, use.first=(param$pet=="first"))
		}

		if (length(cur.reads$pos)) { 
			return(GRanges(cur.chr, IRanges(pmax(1L, cur.reads$pos), 
				pmin(max.len, cur.reads$pos + cur.reads$qwidth - 1L)),
				strand=cur.reads$strand, seqinfo=sqi))
		}
	} else {
		if (param$rescue.pairs) {
			cur.reads <- .rescuePET(bam.file, where=actual.region, dedup=dedup, 
				minq=minq, discard=discard, max.frag=param$max.frag, ext=param$rescue.ext)
		} else {
			cur.reads <- .extractPET(bam.file, where=actual.region, dedup=dedup, 
				minq=minq, discard=discard, max.frag=param$max.frag)
		}

		# Filtering to retain those that don't actually overlap.
		if (length(cur.reads$pos)) { 
			of.interest <- GRanges(cur.chr, IRanges(pmax(1L, cur.reads$pos), 
				pmin(max.len, cur.reads$pos + cur.reads$size - 1L)), seqinfo=sqi)
			keep <- overlapsAny(of.interest, cur.region)
			return(of.interest[keep])
		}
	}
			
	# Returning an empty set, otherwise.
	return(GRanges(seqinfo=sqi))
}

