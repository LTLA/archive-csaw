mergeWindows <- function(regions, tol, sign=NULL, max.width=NULL, ignore.strand=TRUE)
# This function merges the adjacent windows if they lie within 'tol' of each other,
# Any abundance filtering should be done beforehand. Negative values of tol refer
# to a minimum overlap. A value of zero means that the windows must be adjacent
# (i.e. non-overlapping and contiguous).
# 
# written by Aaron Lun
# created 30 July 2013
# last modified 10 February 2015
{
	strs <- strand(regions)
	if (!ignore.strand && length(runValue(strs))!=1) {
		# Strand-specific clustering.
		forward <- as.logical(strs=="+")
		reverse <- as.logical(strs=="-")
		neither <- as.logical(strs=="*")
		ids <- integer(length(regions))
		regs <- NULL
		if (any(forward)) { 
			out <- Recall(regions=regions[forward], tol=tol, sign=sign[forward], max.width=max.width, ignore.strand=TRUE) 
			ids[forward] <- out$id
			strand(out$region) <- "+"
			regs <- out$region
		}
		if (any(reverse)) { 
			out <- Recall(regions=regions[reverse], tol=tol, sign=sign[reverse], max.width=max.width, ignore.strand=TRUE) 
			ids[reverse] <- out$id + length(regs)
			strand(out$region) <- "-"
			regs <- c(regs, out$region)
		}
		if (any(neither)) { 
			out <- Recall(regions=regions[neither], tol=tol, sign=sign[neither], max.width=max.width, ignore.strand=TRUE) 
			ids[neither] <- out$id + length(regs)
			regs <- c(regs, out$region)
		}
		return(list(id=ids, region=regs))
	}

	# Setting up necessary parameters.
	chrs <- as.integer(seqnames(regions))
	starts <- start(regions)
	ends <- end(regions)
	o <- order(chrs, starts, ends)
	if (is.null(sign)) { 
		sign <- logical(length(regions)) 
	} else {
		sign <- sign[o]
	}
	tol <- as.integer(tol)
	max.width <- as.integer(max.width)

	# Running the merge.
	out <- .Call(cxx_merge_windows, chrs[o], starts[o], ends[o], sign, tol, max.width)
	if (is.character(out)) { stop(out) }
	
	# Reporting. Indices correspond with positions in 'clustered'.
	out[[1]][o] <- out[[1]]
	clustered <- GRanges(levels(seqnames(regions))[out[[2]]], IRanges(out[[3]], out[[4]]),
			seqinfo=Seqinfo(seqlevels(regions), seqlengths(regions)))
	return(list(id=out[[1]], region=clustered))
}
