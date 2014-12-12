windowCounts <- function(bam.files, spacing=50, width=spacing, ext=100, shift=0,
	filter=NULL, bin=FALSE, param=readParam())
# Gets counts from BAM files at each position of the sliding window. Applies
# a gentle filter to remove the bulk of window positions with low counts.
# Returns a DGEList with count and total information, as well as a GRanges
# object specifying the window intervals.
# 
# written by Aaron Lun
# ages ago.
# last modified 12 December 2014
{   
	nbam <- length(bam.files)
	paramlist <- .makeParamList(nbam, param)
	extracted.chrs <- .activeChrs(bam.files, paramlist[[1]]$restrict)

	# Processing input parameters.
	if (length(bin)>1L || !is.logical(bin)) { stop("bin must be a logical scalar") }
	if (!bin) { 
		spacing <- as.integer(spacing)
		left <- as.integer(shift)
		right <- as.integer(width) - left - 1L
		ext <- as.integer(ext)
		if (is.null(filter)) { filter <- 5*nbam }
	} else {
		# A convenience flag, which assigns sensible arguments to everything else.
		spacing <- as.integer(width)
		left <- as.integer(shift)
		right <- spacing - 1L - left
		ext <- 1L
		filter <- 1
	}

	# Checking the extension and spacing parameters. We've reparameterised it so
	# that 'left' and 'right' refer to the extension of the window from a nominal
	# 'center' point. This simplifies read counting as we just measure read
	# overlaps to those center points, spaced at regular intervals. 
	if (left >= spacing) { stop("shift must be less than the spacing") }
	if (left < 0L) { stop("shift must be positive") }
	if (left + right < 0L) { stop("width must be a positive integer") }
	if (spacing <= 0L) { stop("spacing must be a positive integer") }
	if (ext <= 0L) { stop("extension width must be a positive integer") }

	# Need to account for the possible loss of a centre point from the front when
	# the shift is non-zero, because the corresponding window is wholly outside the
	# chromosome (i.e., shifted so that the width of the window is before position
	# 1, so width > shift; left+right+1 > left; right+1 > 0; right >= 0).
	at.start <- right >= 0L 
	first.pt <- ifelse(at.start, 1L, spacing+1L)

	# Initializing various collectable containers (non-empty so it'll work if no chromosomes are around).
	totals <- integer(nbam)	
	all.out <- list(matrix(0L, ncol=nbam, nrow=0))
	all.regions <- list(GRanges())
	ix <- 1

	for (i in 1:length(extracted.chrs)) {
		chr <- names(extracted.chrs)[i]
		outlen <- extracted.chrs[i]		
		where <- GRanges(chr, IRanges(1, outlen))

		# Accounting for the possible gain of a centrepoint from the back when
		# shift/left is non-zero, i.e., does the shift bring the next centre point
		# (floor((outlen-1)/spacing)*spacing+1+spacing) under outlen?  [note that 
		# floor - original for a float is equal to the negative remainder].
		at.end <- spacing - (outlen - 1L) %% spacing <= left
		total.pts <- as.integer((outlen-1)/spacing) + at.start + at.end
		outcome <- matrix(0L, total.pts, nbam) 

		for (bf in 1:nbam) {
			curpar <- paramlist[[bf]]
			if (curpar$pet!="both") {
				if (curpar$pet=="none") { 
   					reads <- .extractSET(bam.files[bf], where=where, param=curpar)
				} else {
					reads <- .extractBrokenPET(bam.files[bf], where=where, param=curpar)
				}
				frag.start <- ifelse(reads$strand=="+", reads$pos, reads$pos+reads$qwidth-ext)
				if (length(frag.start)) { frag.start <- pmin(frag.start, outlen) }
				frag.end <- frag.start+ext-1L
			} else {
				if (curpar$rescue.pairs) { 
					out <- .rescuePET(bam.files[bf], where=where, param=curpar)
				} else {
					out <- .extractPET(bam.files[bf], where=where, param=curpar)
				}

				# Only want to record each pair once in a bin, so forcing it to only use the midpoint.
				if (bin) { 
					mid <- as.integer(out$pos + out$size/2)
					frag.end <- frag.start <- mid
				} else { 
					frag.start <- out$pos
					frag.end <- frag.start + out$size - 1L 
				}
			}

			# Extending reads to account for window sizes > 1 bp. The start of each read
			# must be extended by 'right' and the end of each read must be extended by
			# 'left'. We then pull out counts at the specified spacing. We do have to 
			# keep track of whether or not we want to use the first point, though.
			out <- .Call(cxx_get_rle_counts, frag.start-right, frag.end+left, total.pts, spacing, at.start)
			if (is.character(out)) { stop(out) }
			outcome[,bf] <- out
			totals[bf] <- totals[bf]+length(frag.start)
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
		all.regions[[ix]] <- GRanges(chr, IRanges(reg.start, reg.end))
		ix <- ix+1L
	}

	# Generating the remaining GRanges for output (suppressing numerous warnings).
	all.regions <- suppressWarnings(do.call(c, all.regions))
	seqlevels(all.regions) <- names(extracted.chrs)
	seqlengths(all.regions) <- extracted.chrs

	if (is.list(param)) { 
		index <- 1:nbam
	} else {
		index <- 0L
	}
	return(SummarizedExperiment(assays=do.call(rbind, all.out), 
		rowData=all.regions, 
		colData=DataFrame(bam.files=bam.files, totals=totals, param=index),
		exptData=SimpleList(ext=ext, spacing=spacing, width=width, 
			shift=shift, param=param)))
}


