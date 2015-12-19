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
# last modified 19 December 2015
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
            start(where), end(where), param$minq, param$dedup, 
            param$forward, use.first)
    if (is.character(out)) { stop(out) }
    names(out) <- c("forward", "reverse")

	# Filtering by discard regions. Using alignment width so long reads can escape repeats.
    for (i in names(out)) { 
        current <- out[[i]]
        keep <- .discardReads(cur.chr, current[[1]], current[[2]], param$discard)
        current <- lapply(current, "[", keep)
        names(current) <- c("pos", "qwidth", "clip5")
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
    fkeep <- all.sizes$full <= param$max.frag 

    # Reporting output.
    keep <- dkeep & fkeep
    output <- list(pos=left[[1]][keep], size=all.sizes$clipped[keep])
    if (with.reads) {
		output$left <- list(pos=left[[1]][keep], qwidth=left[[2]][keep], strand=rep("+", sum(keep)))
		output$right <- list(pos=right[[1]][keep], qwidth=right[[2]][keep], strand=rep("-", sum(keep)))
	}
	return(output)
}

# Aliases, for convenience.

.getSingleEnd <- .extractSE

.getPairedEnd <- .extractPE

###########################################################

.makeParamList <- function(nbam, param) 
# Converts a readParam object into a list, if it isn't so already.
# 
# written by Aaron Lun
# created 12 December 2014
# last modified 22 July 2015
{
	if (!is.list(param)) { 
		paramlist <- lapply(seq_len(nbam), FUN=function(x) { param })
	} else if (nbam!=length(param)) {
		stop("number of readParam objects is not equal to the number of BAM files")
	} else {
		paramlist <- param
		all.restrict <- paramlist[[1]]$restrict
		for (x in paramlist) {
			if (!identical(x$restrict, all.restrict)) { 
				stop("non-identical restrict settings between readParam objects")
			}		
		}
	}
	return(paramlist)
}	

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

.checkFragments <- function(starts, ends, final, chrlen) 
# Coerces the fragments to the desired 'final.ext', and ensures
# that prior manipulations do not redefine fragment beyond chromosome 
# boundaries (e.g., due to read extension or rescaling).
#
# written by Aaron Lun
# created 13 February 2014
# last modified 14 May 2015
{
	if (!is.na(final)) { 
		remainders <- as.integer((ends - starts + 1L - final)/2)
		if (any(remainders!=0L)) { 
			starts <- starts + remainders
			ends <- ends - remainders
		} 
	}
	if (length(starts)) { starts <- pmin(starts, chrlen) } 
	if (length(ends)) { ends <- pmax(1L, ends) }
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

.extendSE <- function(reads, ext, final, chrlen, retain.strand=FALSE)
# This decides how long to extend reads. The addition of the remainder kicks
# out (or truncates) the fragments to reach the desired 'final.ext'. If 'ext'
# is NA, the read length is used instead.
#
# written by Aaron Lun
# created 12 December 2014
# last modified 19 December 2015
{
	if (is.na(ext)) {   
		frag.start <- c(reads$forward$pos, reads$reverse$pos)
		frag.end <- frag.start + c(reads$forward$qwidth, reads$reverse$qwidth) - 1L	
	} else {
		f.start <- reads$forward$pos
        f.end <- f.start + ext - 1L - reads$forward$clip5 # Cut 5' soft clips out, avoid overextension.
        r.end <- reads$reverse$pos + reads$reverse$qwidth - 1L
        r.start <- r.end - ext + reads$reverse$clip5 + 1L
        frag.start <- c(f.start, r.start)
        frag.end <- c(f.end, r.end)
	}
	out <- .checkFragments(frag.start, frag.end, final=final, chrlen=chrlen)
    if (retain.strand) { out$strand <- Rle(c("+", "-"), c(length(reads$forward[[1]]), length(reads$reverse[[1]]))) }
    return(out)
}

############################################################

.decideStrand <- function(paramlist) 
# Decides what strand we should assign to the output GRanges in the
# SummarizedExperiment object, after counting.
#
# written by Aaron Lun
# created 10 February 2015
{
	getfs <- sapply(paramlist, FUN=function(x) { x$forward })
	if (length(unique(getfs))!=1) {
		warning("unstranded regions used for counts from multiple strands")
		return("*")
	} else if (length(getfs[[1]])==0L) { # Need '[[', if NULL.
		stop("unspecified strandedness")
	}
	if (is.na(getfs[1])) { return("*") }
	else if (getfs[1]) { return("+") }
	else { return("-") }
}

.runningWM <- function(store, x)
# Computes the weighted mean.
{
    if (!length(store)) { store <- c(0, 0) }
    store[1] <- stats::weighted.mean(c(mean(x), store[1]), c(length(x), store[2]))
    store[2] <- store[2] + length(x)
    return(store)
}

.formatColData <- function(bam.files, totals, ext.data, all.pe, all.rlen, paramlist) {
    nbam <- length(bam.files)
    store.ext <- ext.data$ext
    store.rlen <- rep(NA_integer_, nbam)
    for (bf in seq_len(nbam)) {
        if (paramlist[[bf]]$pe=="both") { store.ext[bf] <- as.integer(round(all.pe[[bf]][[1]])) }
        else { store.rlen[bf] <- as.integer(round(all.rlen[[bf]][[1]])) }
    }
    dim(paramlist) <- c(nbam, 1)
    colnames(paramlist) <- "param"
    DataFrame(bam.files=bam.files, totals=totals, ext=store.ext, rlen=store.rlen, paramlist)
}

