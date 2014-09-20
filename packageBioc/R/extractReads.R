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

	# Expanding the extracted region so that it will pull down relevant fragments.
	expand <- 0L
	if (param$pet=="both") {
		expand <- param$max.frag
		if (param$rescue.pairs && expand <= param$rescue.ext) { expand <- param$rescue.ext }
	}
	actual.region <- GRanges(cur.chr, IRanges(max(1L, start(cur.region)-expand),
		min(max.len, end(cur.region)+expand)))

	# Dropping additional reads if required.
	minq <- param$minq
	dedup <- param$dedup
    if (length(param$discard)) { discard <- ranges(param$discard[overlapsAny(param$discard, actual.region)]) }
	else { discard <- NULL }

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
			return(GRanges(cur.chr, IRanges(cur.reads$pos, 
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
			of.interest <- GRanges(cur.chr, IRanges(cur.reads$pos, 
				pmin(max.len, cur.reads$pos + cur.reads$size - 1L)), seqinfo=sqi)
			keep <- overlapsAny(of.interest, cur.region)
			return(of.interest[keep])
		}
	}
			
	# Returning an empty set, otherwise.
	return(GRanges(seqinfo=sqi))
}

