.extractSET <- function(bam, where, param, extras="strand", ...) 
# Extracts single-end read data from a BAM file with removal of unmapped,
# duplicate and poorly mapped/non-unique reads. We also discard reads in the
# specified discard regions. In such cases, the offending reads must be wholly
# within the repeat region.  We use the real alignment width, just in case we
# have very long reads in the alignment that are heavily soft-clipped (i.e., they
# should be reported as within but the read length will put them out).
#
# written by Aaron Lun
# created 8 December 2013
# last modified 12 December 2014
{
	all.fields <- c("pos", "qwidth", extras)
	if (!is.na(param$minq)) { all.fields <- c(all.fields, "mapq") }
	if (length(param$discard)) { all.fields <- c(all.fields, "cigar") }	
	all.fields <- unique(all.fields)
	reads <- scanBam(bam, param=ScanBamParam(what=all.fields,
		which=where, flag=scanBamFlag(isUnmappedQuery=FALSE, 
		isDuplicate=ifelse(param$dedup, FALSE, NA), ...)))[[1]]
   
	# Filtering by MAPQ.
	if (!is.na(param$minq)) { 
		keep <- reads$mapq >= param$minq & !is.na(reads$mapq) 
		if (!"mapq" %in% extras) { reads$mapq <- NULL }
		for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
	}
	
	# Filtering by discard regions. Using alignment width so long reads can escape repeats.
	if (length(param$discard)) {
		relevant <- seqnames(param$discard)==as.character(seqnames(where[1]))
		if (any(relevant)) { 
			awidth <- cigarWidthAlongReferenceSpace(reads$cigar)
			keep <- !overlapsAny(IRanges(reads$pos, reads$pos+awidth-1L), 
				ranges(param$discard[relevant]), type="within")
			for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
		}
		if (!"cigar" %in% extras) { reads$cigar <- NULL }
	}
	return(reads)
}

.extractPET <- function(bam.file, where, param)
# A function to extract PET data for a particular chromosome. Synchronisation
# is expected.  We avoid sorting by name  as it'd mean we have to process the
# entire genome at once (can't go chromosome-by-chromosome).  This probably
# results in increased memory usage across the board, and it doesn't fit in
# well with the rest of the pipelines which assume coordinate sorting.
# 
# written by Aaron Lun
# created 8 December 2013
# last modified 12 December 2012
{
	reads <- .extractSET(bam.file, extras=c("qname", "flag"), where=where, 	
		param=param, isPaired=TRUE, hasUnmappedMate=FALSE)
	.yieldInterestingBits(reads, max(end(where)), max.frag=param$max.frag)
}

.extractBrokenPET <- function(bam.file, where, param)
# A function to extract PET data, but as single-end data (i.e. only using one
# of the reads).  Useful when paired-end data has gone completely off the
# rails.
{
	use.first <- param$pet=="first"
	.extractSET(bam.file, where=where, param=param,  
		isPaired=TRUE, isFirstMateRead=use.first, isSecondMateRead=!use.first)
}

.rescuePET <- function(bam.file, where, param)
# A function to extract PET data where possible, but to rescue those that
# are invalid by using them as single-end data with read extension. Those
# reads that form invalid pairs are broken up and the read with the better
# MAPQ is chosen. Any single read (due to filtering or whatever) is used as-is.
# Interchromosomal pairs get counted once on each chromosome.
#
# written by Aaron Lun
# created 13 May 2014
# last modified 12 December 2012
{
	reads <- .extractSET(bam.file, extras=c("qname", "flag", "mapq"), where=where, 
		param=param, isPaired=TRUE)
	output <- .yieldInterestingBits(reads, max(end(where)), diag=TRUE, max.frag=param$max.frag)

	# Figuring out which reads do pair up.
	is.first <- bitwAnd(reads$flag, 0x40) != 0L
	nok.first <- !output$is.ok & is.first
	nok.second <- !output$is.ok & !is.first
	corresponding <- match(reads$qname[nok.first], reads$qname[nok.second])
	first.paired <- which(!is.na(corresponding))
	second.paired <- corresponding[first.paired]
	
	# Picking the unique reads, or the better read from the each pair.
	additor <- logical(length(reads$flag))
	additor[nok.first][-first.paired] <- TRUE
	additor[nok.second][-second.paired] <- TRUE
	is.better <- reads$mapq[nok.first][first.paired] > reads$mapq[nok.second][second.paired]
	additor[nok.first][first.paired][is.better] <- TRUE
	additor[nok.second][second.paired][!is.better] <- TRUE

	# Returning the loot.
	return( list( pos=c(output$pos, ifelse(bitwAnd(reads$flag[additor], 0x10)==0L, reads$pos[additor], 
					reads$pos[additor]+reads$qwidth[additor]-param$ext)),
		      size=c(output$size, rep(param$ext, sum(additor))) ) )
}

.makeParamList <- function(nbam, param) 
# Converts a readParam object into a list, if it isn't so already.
# 
# written by Aaron Lun
# created 12 December 2014
{
	if (!is.list(param)) { 
		paramlist <- lapply(1:nbam, FUN=function(x) { param })
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

.extendSE <- function(reads, chrlen, param)
# This decides how long to extend reads. There's two parameters here; how long
# to extend the read, and what the final fragment size should be. These are
# separate parameters, to allow synchronization of different libraries to the
# same fragment length.
#
# written by Aaron Lun
# created 12 December 2014
{
	frag.start <- ifelse(reads$strand=="+", reads$pos, reads$pos + reads$qwidth - param$ext)
	if (length(frag.start)) { frag.start <- pmin(frag.start, chrlen) } # to ensure the read is still counted when ext < qwidth.
	frag.end <- frag.start + param$ext - 1L
	return(list(start=frag.start, end=frag.end))
}
