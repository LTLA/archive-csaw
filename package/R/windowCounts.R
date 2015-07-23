windowCounts <- function(bam.files, spacing=50, width=spacing, ext=100, shift=0,
	filter=10, bin=FALSE, param=readParam())
# Gets counts from BAM files at each position of the sliding window. Applies
# a gentle filter to remove the bulk of window positions with low counts.
# Returns a RangedSummarizedExperiment object with counts and genomic
# coordinates.
# 
# written by Aaron Lun
# created 5 April 2012
# last modified 22 July 2015
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
	} else {
		# A convenience flag, which assigns sensible arguments to everything else.
		spacing <- as.integer(width)
		left <- as.integer(shift)
		right <- spacing - 1L - left
		ext <- 1L
		final.ext <- NA
		filter <- min(1, filter)
	}
	ext.data <- .collateExt(nbam, ext)

	# Checking the extension and spacing parameters. We've reparameterised it so
	# that 'left' and 'right' refer to the extension of the window from a nominal
	# 'center' point. This simplifies read counting as we just measure read
	# overlaps to those center points, spaced at regular intervals. 
	if (left >= spacing) { stop("shift must be less than the spacing") }
	if (left < 0L) { stop("shift must be positive") }
	if (left + right < 0L) { stop("width must be a positive integer") }
	if (spacing <= 0L) { stop("spacing must be a positive integer") }

	# Need to account for the possible loss of a centre point from the front when
	# the shift is non-zero, because the corresponding window is wholly outside the
	# chromosome (i.e., shifted so that the width of the window is before position
	# 1; can't happen if width > shift; left+right+1 > left; right+1 > 0; right >= 0).
	at.start <- right >= 0L 
	first.pt <- ifelse(at.start, 1L, spacing+1L)

	# Initializing various collectable containers (non-empty so it'll work if no chromosomes are around).
	totals <- integer(nbam)	
	all.out <- list(matrix(0L, ncol=nbam, nrow=0))
	all.regions <- list(GRanges())
	ix <- 1

	for (i in seq_along(extracted.chrs)) {
		chr <- names(extracted.chrs)[i]
		outlen <- extracted.chrs[i]		
		where <- GRanges(chr, IRanges(1, outlen))

		# Accounting for the possible gain of a centrepoint from the back when
		# shift/left is non-zero, i.e., does the shift bring the next centre point
		# (floor((outlen-1)/spacing)*spacing+1+spacing) under outlen?  [note that 
		# floor - original for (outlen-1) is equal to the negative remainder].
		at.end <- spacing - (outlen - 1L) %% spacing <= left
		total.pts <- as.integer((outlen-1)/spacing) + at.start + at.end
		outcome <- matrix(0L, total.pts, nbam) 

		for (bf in seq_len(nbam)) {
			curpar <- paramlist[[bf]]
			if (curpar$pe!="both") {
   				reads <- .getSingleEnd(bam.files[bf], where=where, param=curpar)
				extended <- .extendSE(reads, ext=ext.data$ext[bf], final=ext.data$final, chrlen=outlen)
				frag.start <- extended$start
				frag.end <- extended$end
			} else {
				out <- .getPairedEnd(bam.files[bf], where=where, param=curpar)
				if (bin) { 
					# Only want to record each pair once in a bin, so forcing it to only use the midpoint.
					mid <- as.integer(out$pos + out$size/2)
					frag.end <- frag.start <- mid
				} else {
					checked <- .checkFragments(out$pos, out$pos+out$size-1L, final=ext.data$final, chrlen=outlen)
					frag.start <- checked$start
					frag.end <- checked$end
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
	strand(all.regions) <- .decideStrand(paramlist)

	dim(paramlist) <- c(nbam, 1)
	colnames(paramlist) <- "param"
	return(SummarizedExperiment(assays=SimpleList(counts=do.call(rbind, all.out)), 
		rowRanges=all.regions, 
		colData=DataFrame(bam.files=bam.files, totals=totals, ext=ext.data$ext, paramlist),
		metadata=list(spacing=spacing, width=width, shift=shift, 
			final.ext=ifelse(bin, 1L, ext.data$final)))) # For getWidths with paired-end binning.
}

