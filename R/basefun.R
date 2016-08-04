.extractSE <- function(bam.file, where, param) 
# Extracts single-end read data from a BAM file with removal of unmapped,
# duplicate and poorly mapped/non-unique reads. We also discard reads in the
# specified discard regions. In such cases, the offending reads must be wholly
# within the repeat region.  We use the real alignment width, just in case we
# have very long reads in the alignment that are heavily soft-clipped (i.e., they
# should be reported as within but the read length will put them out).
#
# written by Aaron Lun
# created 8 December 2013
# last modified 20 December 2015
{
    cur.chr <- as.character(seqnames(where)) 
    bam.file <- path.expand(bam.file)
    bam.index <- paste0(bam.file, ".bai")

	if (length(param$forward)==0L) { stop("read strand extraction must be specified") }
    if (param$pe=="first") {
        use.first <- TRUE
    } else if (param$pe=="second") { 
        use.first <- FALSE
    } else {
        use.first <- NA
    }

    out <- .Call(cxx_extract_single_data, bam.file, bam.index, cur.chr,
            start(where), end(where), param$minq, param$dedup, param$forward, use.first)
    if (is.character(out)) { stop(out) }
    
    names(out) <- c("forward", "reverse")
    names(out$forward) <- names(out$reverse) <- c("pos", "qwidth")

	# Filtering by discard regions. Using alignment width so long reads can escape repeats.
    for (i in names(out)) { 
        current <- out[[i]]
        keep <- .discardReads(cur.chr, current[[1]], current[[2]], param$discard)
        current <- lapply(current, "[", keep)
        out[[i]] <- current
    }
    return(out)
}

.discardReads <- function(chr, pos, alen, discard) {
    if (!length(pos)) { 
        return(logical(0)) 
    }
    relevant <- seqnames(discard)==chr
    if (any(relevant)) { 
		keep <- !overlapsAny(IRanges(pos, pos+alen-1L), ranges(discard)[relevant], type="within")
    } else {
        keep <- !logical(length(pos))
    }
    return(keep)
}

.extractPE <- function(bam.file, where, param, with.reads=FALSE, diagnostics=FALSE)
# A function to extract PE data for a particular chromosome. Synchronisation
# is expected.  We avoid sorting by name  as it'd mean we have to process the
# entire genome at once (can't go chromosome-by-chromosome).  This probably
# results in increased memory usage across the board, and it doesn't fit in
# well with the rest of the pipelines which assume coordinate sorting.
# 
# written by Aaron Lun
# created 8 December 2013
# last modified 16 December 2015
{
    cur.chr <- as.character(seqnames(where)) 
    bam.file <- path.expand(bam.file)
    bam.index <- paste0(bam.file, ".bai")
    out <- .Call(cxx_extract_pair_data, bam.file, bam.index, cur.chr,
            start(where), end(where), param$minq, param$dedup, diagnostics)
    if (is.character(out)) { stop(out) }

    if (diagnostics) {
        names(out) <- c("forward", "reverse", "total", "single", "ufirst", "usecond", "one.mapped", "ifirst", "isecond")
        return(out)
    }
    left <- out[[1]]
    right <- out[[2]]

    # Filtering by discard.
    dlkeep <- .discardReads(cur.chr, left[[1]], left[[2]], param$discard)
    drkeep <- .discardReads(cur.chr, right[[1]], right[[2]], param$discard)
    dkeep <- dlkeep & drkeep

    # Filtering by maximum fragment size.
    all.sizes <- .getFragmentSizes(left, right)
    fkeep <- all.sizes <= param$max.frag 

    # Reporting output.
    keep <- dkeep & fkeep
    output <- list(pos=left[[1]][keep], size=all.sizes[keep])
    if (with.reads) {
		output$forward <- lapply(left, "[", keep)
		output$reverse <- lapply(right, "[", keep)
        names(output$forward) <- names(output$reverse) <- c("pos", "qwidth")
    }
	return(output)
}

# Aliases, for convenience.

.getSingleEnd <- .extractSE

.getPairedEnd <- .extractPE

###########################################################

.activeChrs <- function(bam.files, restrict) 
# Processes the incoming data; checks that bam headers are all correct,
# truncates the list according to 'restrict'.
# 
# written by Aaron Lun
# created 12 December 2014
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
    if (length(restrict)) { originals <- originals[names(originals) %in% restrict] }
	return(originals)
}

.boundIntervals <- function(x, chrlen) {
    x[x<1L] <- 1L
    x[x>chrlen] <- chrlen
    return(x)
}

.coerceFragments <- function(starts, ends, final, chrlen) 
# Coerces the fragments to the desired 'final.ext', and ensures
# that prior manipulations do not redefine fragment beyond chromosome 
# boundaries (e.g., due to read extension or rescaling).
#
# written by Aaron Lun
# created 13 February 2014
# last modified 21 December 2015
{
	if (!is.na(final)) { 
		remainders <- as.integer((ends - starts + 1L - final)/2)
		if (any(remainders!=0L)) { 
			starts <- starts + remainders
			ends <- ends - remainders
		} 
	}
    starts <- .boundIntervals(starts, chrlen)
    ends <- .boundIntervals(ends, chrlen)
	return(list(start=starts, end=ends)) 
}

############################################################

.collateExt <- function(nbam, ext)
# Collates the extension parameters into a set of ext and remainder values.
# The idea is to extend each read directionally to 'ext', and then extend in
# both directions by 'remainder' to reach the desired fragment length.
# 
# written by Aaron Lun
# created 12 December 2014
# last modified 16 December 2015
{
    if (!is.vector(ext)) {
        if (length(ext)!=2L) {
            stop("'ext' must be a list of length 2")
        } else {
            final.ext <- unique(as.integer(round(ext[[2]])))
        	if (length(final.ext)!=1L || (!is.na(final.ext) && final.ext <= 0L)) { 
                stop("final extension length must be a positive integer or NA") 
            }
        }
        ext <- ext[[1]]
    } else {
        final.ext <- NA_integer_
    }
	
	if (length(ext)==1L) { 
		ext <- rep(ext, nbam)
	} else if (length(ext)!=nbam) {
		stop("length of extension vector is not consistent with number of libraries")
	}
	ext <- as.integer(round(ext))
	if (any(!is.na(ext) & ext <= 0L)) { stop("extension length must be NA or a positive integer") }

	list(ext=ext, final=final.ext)
}

.extendSEdir <- function(reads, ext, final, chrlen, forward=TRUE) {
    if (is.na(ext)) { 
        start <- reads$pos
        end <- reads$pos + reads$qwidth -1L
    } else {
        if (forward) {
            start <- reads$pos
            end <- reads$pos + ext - 1L
        } else {
            end <- reads$pos + reads$qwidth - 1L
            start <- end - ext + 1L
        }
    }
	out <- .coerceFragments(start, end, final=final, chrlen=chrlen)
    return(out)
}

.extendSE <- function(reads, ext, final, chrlen)
# This decides how long to extend reads. The addition of the remainder kicks
# out (or truncates) the fragments to reach the desired 'final.ext'. If 'ext'
# is NA, the read length is used instead.
#
# written by Aaron Lun
# created 12 December 2014
# last modified 21 December 2015
{
    fout <- .extendSEdir(reads$forward, ext, final, chrlen, forward=TRUE)
    rout <- .extendSEdir(reads$reverse, ext, final, chrlen, forward=FALSE)
    mapply(c, fout, rout, SIMPLIFY=FALSE)
}

############################################################

.decideStrand <- function(param) 
# Decides what strand we should assign to the output GRanges in the
# SummarizedExperiment object, after counting.
{
	getfs <- param$forward
    if (length(getfs)==0L) { 
		stop("unspecified strandedness")
	} else if (is.na(getfs)) { 
        return("*") 
    } else if (getfs) { 
        return("+") 
    } else { 
        return("-") 
    }
}

.runningWM <- function(store, x)
# Computes the weighted mean.
{
    if (!length(store)) { store <- c(0, 0) }
    store[1] <- weighted.mean(c(mean(x), store[1]), c(length(x), store[2]))
    store[2] <- store[2] + length(x)
    return(store)
}

.formatColData <- function(bam.files, totals, ext.data, all.pe, all.rlen, param) {
    nbam <- length(bam.files)
    store.ext <- ext.data$ext
    store.rlen <- rep(NA_integer_, nbam)
    for (bf in seq_len(nbam)) {
        if (param$pe=="both") { store.ext[bf] <- as.integer(round(all.pe[[bf]][[1]])) }
        else { store.rlen[bf] <- as.integer(round(all.rlen[[bf]][[1]])) }
    }
    DataFrame(bam.files=bam.files, totals=totals, ext=store.ext, rlen=store.rlen)
}

############################################################

.toGRanges <- function(x) 
# Converts the input to a GRanges, if it wasn't before.
{
    if (is(x, "RangedSummarizedExperiment")) { 
        x <- rowRanges(x)
    } else if (!is(x, "GenomicRanges")) {
        stop("'x' must be a RangedSummarizedExperiment or GRanges object")
    }
    return(x)
}
