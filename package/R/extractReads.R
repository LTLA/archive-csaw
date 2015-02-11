extractReads <- function(cur.region, bam.file, param=readParam())
# Exactly as specified. Takes a region and plots it in bimodal mode, with
# options for duplicate removal, mapping quality enhancement, colorization,
# and PE manipulation.
#
# written by Aaron Lun
# created 1 September 2014
# last modified 10 February 2015
{
    if (length(cur.region)!=1L) { stop("exactly one range is required for plotting") }
	if (as.logical(strand(cur.region)!="*")) { warning("strandedness of region will be ignored, use param$forward instead") }

    chrs <- scanBamHeader(bam.file)[[1]][[1]]
	cur.chr <- as.character(seqnames(cur.region)[1])
	if (length(param$restrict) && ! cur.chr %in% param$restrict) { stop("current chromosome not in restricted subset") }
	if (! cur.chr %in% names(chrs)) { stop("cannot find current chromosome in the BAM file header") }
	max.len <- chrs[[cur.chr]]
	sqi <- Seqinfo(cur.chr, max.len)

	# Extracting all-of-chromosome for paired-end rescue, as you need to find the read with the higher MAPQ.
	expand <- 0L
	if (param$pe=="both" && !is.na(param$rescue.ext)) {
		actual.region <- GRanges(cur.chr, IRanges(1L, max.len)) 
	} else {
		if (param$pe=="both") {
			actual.region <- GRanges(cur.chr, IRanges(max(1L, start(cur.region)-param$max.frag),
				min(max.len, end(cur.region)+param$max.frag)))
		} else {
			actual.region <- cur.region
			seqlevels(actual.region) <- cur.chr
		}
	}

	# Pulling out reads from a region and setting up coverage RLE's.
	if (param$pe!="both") {
		if (param$pe=="none") { 
			cur.reads <- .extractSE(bam.file, where=actual.region, param=param)
		} else {
			cur.reads <- .extractBrokenPE(bam.file, where=actual.region, param=param)
		}

		if (length(cur.reads$pos)) { 
			return(GRanges(cur.chr, IRanges(pmax(1L, cur.reads$pos), 
				pmin(max.len, cur.reads$pos + cur.reads$qwidth - 1L)),
				strand=cur.reads$strand, seqinfo=sqi))
		}
	} else {
		if (!is.na(param$rescue.ext)) {
			cur.reads <- .rescuePE(bam.file, where=actual.region, param=param)
		} else {
			cur.reads <- .extractPE(bam.file, where=actual.region, param=param)
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

