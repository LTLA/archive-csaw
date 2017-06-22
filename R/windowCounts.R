windowCounts <- function(bam.files, spacing=50, width=spacing, ext=100, 
    shift=0, filter=10, bin=FALSE, param=readParam())
# Gets counts from BAM files at each position of the sliding window. Applies
# a gentle filter to remove the bulk of window positions with low counts.
# Returns a RangedSummarizedExperiment object with counts and genomic
# coordinates.
# 
# written by Aaron Lun
# created 5 April 2012
# last modified 22 June 2017
{   
	nbam <- length(bam.files)
	extracted.chrs <- .activeChrs(bam.files, param$restrict)

	# Processing input parameters.
    shift <- as.integer(shift)
    width <- as.integer(width)
    spacing <- as.integer(spacing)
	if (length(bin)>1L || !is.logical(bin)) { 
        stop("bin must be a logical scalar") 
    }
	if (bin) { # A convenience flag for binning.
		spacing <- width
		ext <- 1L
		filter <- min(1, filter)
	}
	ext.data <- .collateExt(nbam, ext)

	# Checking the values of the parameters. 
	if (shift >= spacing) { stop("shift must be less than the spacing") }
	if (shift < 0L) { stop("shift must be positive") }
	if (width <= 0L) { stop("width must be a positive integer") }
	if (spacing <= 0L) { stop("spacing must be a positive integer") }

	# Initializing various collectable containers (non-empty so it'll work if no chromosomes are around).
	totals <- integer(nbam)
    nchrs <- length(extracted.chrs)    
    all.out <- rep(list(matrix(0L, ncol=nbam, nrow=0)), nchrs) # ensure that the class is right, even if nothing is added.
	all.regions <- rep(list(GRanges()), nchrs)
    all.lengths <- rep(list(vector("list", nchrs)), nbam)

	for (i in seq_len(nchrs)) { 
		chr <- names(extracted.chrs)[i]
		outlen <- extracted.chrs[i]		
		where <- GRanges(chr, IRanges(1, outlen))

        # Need to account for the possible loss of a window from the front when
        # the shift is non-zero, because the corresponding window is wholly outside the
        # chromosome (i.e., shifted so that the width of the window is before position 1).
        at.start <- shift < width 

		# Accounting for the possible gain of a centrepoint from the back when
		# shift is non-zero, i.e., does the shift bring the next window start under outlen?
        # i.e., 1L - shift + spacing * X <= outlen, solve for the largest integer X.
        total.pts <- as.integer(floor((outlen + shift - 1L)/spacing))
        total.pts <- total.pts + at.start # if the extra point exists at the start.
		
        # Parallelized loading.
        bp.out <- bpmapply(FUN=.window_counts, bam.file=bam.files, init.ext=ext.data$ext, 
                           MoreArgs=list(where=where, param=param, 
                                         final.ext=ext.data$final, outlen=outlen, bin=bin, 
                                         shift=shift, width=width, spacing=spacing, 
                                         total.pts=total.pts, at.start=at.start),
                           BPPARAM=param$BPPARAM, SIMPLIFY=FALSE)
        
        outcome <- matrix(0L, total.pts, nbam)
        for (bf in seq_along(bp.out)) {
            outcome[,bf] <- bp.out[[bf]]$counts
            totals[bf] <- totals[bf] + bp.out[[bf]]$totals
            all.lengths[[bf]][[i]] <- bp.out[[bf]]$lengths
        }

		# Filtering on row sums (for memory efficiency). Note that redundant windows
		# are avoided by enforcing 'shift < spacing', as a window that is shifted to 
		# cover the whole chromosome cannot be spaced to a new position where it still 
		# covers the whole chromosome. The next one must be inside the chromosome.
		keep <- rowSums(outcome)>=filter 
		if (!any(keep)) { 
            next 
        } else if (!all(keep)) { 
            outcome <- outcome[keep,,drop=FALSE] 
        }
		all.out[[i]] <- outcome

        # Defining genomic coordinates.
		center <- (which(keep) - 1L) * spacing + 1L + ifelse(at.start, 0L, spacing) 
		reg.start <- center - shift
		reg.end <- pmin(outlen, reg.start + width - 1L)
        reg.start <- pmax(1L, reg.start)
		all.regions[[i]] <- GRanges(chr, IRanges(reg.start, reg.end))
	}

	# Generating the remaining GRanges for output (suppressing numerous warnings).
	all.regions <- suppressWarnings(do.call(c, all.regions))
	seqlevels(all.regions) <- names(extracted.chrs)
	seqlengths(all.regions) <- extracted.chrs
	strand(all.regions) <- .decideStrand(param)

	return(SummarizedExperiment(assays=SimpleList(counts=do.call(rbind, all.out)), 
		rowRanges=all.regions, 
		colData=.formatColData(bam.files, totals, ext.data, all.lengths, param),
		metadata=list(spacing=spacing, width=width, shift=shift, bin=bin, 
            param=param, final.ext=ifelse(bin, 1L, ext.data$final)))) # For getWidths with paired-end binning.
}

.window_counts <- function(bam.file, where, param, 
                           init.ext, final.ext, outlen, bin, 
                           shift, width, spacing, 
                           total.pts, at.start) {

    if (param$pe!="both") {
        reads <- .getSingleEnd(bam.file, where=where, param=param)
        extended <- .extendSE(reads, ext=init.ext, final=final.ext, chrlen=outlen)
        frag.start <- extended$start
        frag.end <- extended$end

        rlengths <- cbind(c(mean(reads$forward$qwidth), mean(reads$reverse$qwidth)),
                       c(length(reads$forward$qwidth), length(reads$reverse$qwidth)))
    } else {
        out <- .getPairedEnd(bam.file, where=where, param=param)
        if (bin) { 
            # Only want to record each pair once in a bin, so forcing it to only use the midpoint.
            mid <- as.integer(out$pos + out$size/2)
            frag.end <- frag.start <- mid
        } else {
            checked <- .coerceFragments(out$pos, out$pos+out$size-1L, final=final.ext, chrlen=outlen)
            frag.start <- checked$start
            frag.end <- checked$end
        }
        rlengths <- c(mean(out$size), length(out$size))
    }

    # Windows are parametrized in terms of extension from the window start.
    # Non-unity window width corresponds to extension to the right.
    # Shifting causes it to be extended to the left, at the expense of the right.
    right <- width - shift - 1L
    left <- shift

    # We make it even simpler by extending the reads, which means that we don't 
    # have to actually change the window starts. The *start* of each read must be 
    # extended by 'right' and the *end* of each read must be extended by 'left'. 
    # (This mirrors the extension of the windows.)
    frag.start <- frag.start - right
    frag.end <- frag.end + left    
    
    # We pull out counts at the specified spacing. We do have to keep track of 
    # whether or not we want to use the first point, though.
    out <- .Call(cxx_get_rle_counts, frag.start, frag.end, total.pts, spacing, at.start)
    if (is.character(out)) { stop(out) }
    return(list(counts=out, totals=length(frag.start), lengths=rlengths))
}

