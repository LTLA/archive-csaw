.extractSE <- function(bam, where, param, extras="strand", discard=NULL, ...) 
# Extracts single-end read data from a BAM file with removal of unmapped,
# duplicate and poorly mapped/non-unique reads. We also discard reads in the
# specified discard regions. In such cases, the offending reads must be wholly
# within the repeat region.  We use the real alignment width, just in case we
# have very long reads in the alignment that are heavily soft-clipped (i.e., they
# should be reported as within but the read length will put them out).
#
# written by Aaron Lun
# created 8 December 2013
# last modified 2 July 2015
{
	all.fields <- c("pos", "qwidth", extras)
	if (length(discard)) { all.fields <- c(all.fields, "cigar") }	
	all.fields <- unique(all.fields)

	if (length(param$forward)==0L) { stop("read strand extraction must be specified") }
	if (length(where)!=1L) { stop("extraction not supported for multiple ranges at once") }
	reads <- scanBam(bam, param=ScanBamParam(what=all.fields, which=where, mapqFilter=param$minq,
		flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=ifelse(param$dedup, FALSE, NA), 
		isMinusStrand=ifelse(is.na(param$forward), NA, !param$forward), ...)))[[1]]
  	
	# Filtering by discard regions. Using alignment width so long reads can escape repeats.
	if (length(discard)) {
	    awidth <- cigarWidthAlongReferenceSpace(reads$cigar)
        keep <- .discardReads(as.character(seqnames(where)), reads$pos, awidth, param$discard)
		if (!"cigar" %in% extras) { reads$cigar <- NULL }
		for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
	}
	return(reads)
}

.discardReads <- function(chr, pos, alen, discard) {
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
# last modified 15 December 2015
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
    if (!is.na(param$max.frag)) { 
        fkeep <- all.sizes$full <= param$max.frag 
    } else {
        fkeep <- !logical(length(all.sizes$full))
    }

    # Reporting output.
    keep <- dkeep & fkeep
    output <- list(pos=left[keep], size=all.sizes$clipped[keep])
    if (with.reads) {
		left <- list(pos=left[[1]][keep], qwidth=left[[2]][keep], strand=rep("+", sum(keep)))
		right <- list(pos=right[[1]][keep], qwidth=right[[2]][keep], strand=rep("-", sum(keep)))
		output$left <- left
 	   	output$right <- right
	}
	return(output)
}

.extractBrokenPE <- function(bam.file, where, param)
# A function to extract PE data, but as single-end data (i.e. only using one
# of the reads). Useful when paired-end data has gone completely off the rails.
{
	use.first <- param$pe=="first"
	.extractSE(bam.file, where=where, param=param,  
		isPaired=TRUE, isFirstMateRead=use.first, isSecondMateRead=!use.first)
}

.getSingleEnd <- function(bam, where, param) {
	if (param$pe=="none") { 
		reads <- .extractSE(bam, where=where, param=param)
    } else {
		reads <- .extractBrokenPE(bam, where=where, param=param)
	}
	return(reads)
}	

.getPairedEnd <- function(bam, where, param, with.reads=FALSE) {
	.extractPE(bam, where=where, param=param, with.reads=with.reads)
}

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

.collateExt <- function(nbam, ext)
# Collates the extension parameters into a set of ext and remainder values.
# The idea is to extend each read directionally to 'ext', and then extend in
# both directions by 'remainder' to reach the desired fragment length.
# 
# written by Aaron Lun
# created 12 December 2014
# last modified 13 February 2015
{
	final.ext <- attributes(ext)$final.ext # Do this, before attributes are lost.
	if (is.null(final.ext)) { final.ext <- NA }
	final.ext <- as.integer(final.ext)
	if (length(final.ext)!=1L || (!is.na(final.ext) && final.ext <= 0L)) { 
		stop("final extension length must be a positive integer or NA") 
	}
	
	if (length(ext)==1L) { 
		ext <- rep(ext, nbam)
	} else if (length(ext)!=nbam) {
		stop("length of extension vector is not consistent with number of libraries")
	}
	ext <- as.integer(ext)
	if (any(!is.na(ext) & ext <= 0L)) { stop("extension length must be NA or a positive integer") }

	list(ext=ext, final=final.ext)
}

.extendSE <- function(reads, ext, final, chrlen)
# This decides how long to extend reads. The addition of the remainder kicks
# out (or truncates) the fragments to reach the desired 'final.ext'. If 'ext'
# is NA, the read length is used instead.
#
# written by Aaron Lun
# created 12 December 2014
# last modified 14 May 2015
{
	if (is.na(ext)) {   
		frag.start <- reads$pos
		frag.end <- frag.start + reads$qwidth - 1L	
	} else {
		is.reverse <- reads$strand!="+"
		frag.start <- reads$pos
		if (any(is.reverse)) { 
			frag.start[is.reverse] <- reads$pos[is.reverse] + reads$qwidth[is.reverse] - ext
		}
		frag.end <- frag.start + ext - 1L
	}
	.checkFragments(frag.start, frag.end, final=final, chrlen=chrlen)
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

