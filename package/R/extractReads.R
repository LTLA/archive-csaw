extractReads <- function(bam.file, region, ext=NA, param=readParam(), as.reads=FALSE)
# Exactly as specified. Takes a region and plots it in bimodal mode, with
# options for duplicate removal, mapping quality enhancement, colorization,
# and PE manipulation.
#
# written by Aaron Lun
# created 1 September 2014
# last modified 22 July 2015
{
    if (length(region)!=1L) { stop("exactly one range is required for read extraction") }
	if (as.logical(strand(region)!="*")) { warning("strandedness of region will be ignored, use param$forward instead") }
	ext.data <- .collateExt(1, ext)

    chrs <- scanBamHeader(bam.file)[[1]][[1]]
	cur.chr <- as.character(seqnames(region)[1])
	if (length(param$restrict) && ! cur.chr %in% param$restrict) { warning("current chromosome not in restricted subset") }
	if (! cur.chr %in% names(chrs)) { stop("cannot find current chromosome in the BAM file header") }
	max.len <- chrs[[cur.chr]]
	sqi <- Seqinfo(cur.chr, max.len)

	# Extracting all-of-chromosome for paired-end rescue, as you need to find the read with the higher MAPQ.
	# Otherwise, kicking out the region by 'max.frag' or 'ext' to guarantee capture of all participants.
	if (param$pe=="both" && !param$fast.pe && .needsRescue(param)) { 
		actual.region <- GRanges(cur.chr, IRanges(1L, max.len)) 
	} else {
		max.ext <- suppressWarnings(max(ext.data$ext, ext.data$final, param$max.frag, na.rm=TRUE))
		if (max.ext < 0L) { max.ext <- 0L }
		actual.region <- GRanges(cur.chr, IRanges(max(1L, start(region)-max.ext),
			min(max.len, end(region)+max.ext)))
	}

	# Pulling out reads from a region and setting up coverage RLE's.
	if (param$pe!="both") {
		read.data <- .getSingleEnd(bam.file, where=actual.region, param=param)
		stranded <- read.data$strand

		if (length(stranded)) { 
			cur.reads <- .extendSE(read.data, ext=ext.data$ext[1], final=ext.data$final, chrlen=max.len)
			of.interest <- GRanges(cur.chr, IRanges(pmax(1L, cur.reads$start), pmin(max.len, cur.reads$end)), strand=stranded, seqinfo=sqi)

			if (max.ext) { 
 			   	# Filtering to retain those extended reads that actually overlap.	
				keep <- overlapsAny(of.interest, region)
				of.interest <- of.interest[keep]
			}
			return(of.interest)
		}
	} else {
		frag.data <- .getPairedEnd(bam.file, where=actual.region, param=param, with.reads=as.reads)
		cur.frags <- .checkFragments(frag.data$pos, frag.data$pos + frag.data$size - 1L, final=ext.data$final, chrlen=max.len)

		if (length(cur.frags$start)) { 
			# Filtering to retain those fragments that actually overlap.
			of.interest <- GRanges(cur.chr, IRanges(pmax(1L, cur.frags$start), pmin(max.len, cur.frags$end)), seqinfo=sqi)
			keep <- overlapsAny(of.interest, region)
			if (!as.reads) { return(of.interest[keep]) }

			# Reporting the individual reads, if requested.
			npairs <- length(frag.data$left$pos)
			if (npairs) { 
				left <- suppressWarnings(GRanges(cur.chr, IRanges(frag.data$left$pos, frag.data$left$pos+frag.data$left$qwidth-1L), 
					seqinfo=sqi, strand=frag.data$left$strand))
				right <- suppressWarnings(GRanges(cur.chr, IRanges(frag.data$right$pos, frag.data$right$pos+frag.data$right$qwidth-1L),
					seqinfo=sqi, strand=frag.data$right$strand))

				pairdex <- seq_len(npairs) # first lot of fragments correspond to proper pairs.
				left <- left[keep[pairdex]]
				right <- right[keep[pairdex]]
				left$pair <- right$pair <- seq_along(left)
				reads <- suppressWarnings(c(left, right))
			} else { reads <- NULL }

			nrescue <- length(frag.data$rescued$pos)
			if (nrescue) { 
				rescued <- suppressWarnings(GRanges(cur.chr, IRanges(frag.data$rescued$pos, frag.data$rescued$pos+frag.data$rescued$qwidth), 
					seqinfo=sqi, strand=frag.data$rescued$strand, pair=integer(length(frag.data$rescued$pos))))
				rescued <- rescued[keep[npairs + seq_len(nrescue)]] # second lot correspond to rescued reads.
				reads <- suppressWarnings(c(reads, rescued))
			}

			return(reads)
		}
	}
			
	# Returning an empty set, otherwise.
	return(GRanges(seqinfo=sqi))
}

