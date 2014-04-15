getPETSizes <- function(bam.file, dedup=FALSE, minq=0, restrict=NULL, discard=NULL) 
# This function takes a BAM file and reads it to parse the size of the PET fragments. It then
# returns a vector of sizes which can be plotted for diagnostics. The length of the vector
# will also tell you how many read pairs were considered valid. The total number of reads, the
# number of singletons and the number of interchromosomal pairs is also reported.
# 
# written by Aaron Lun
# a long long time ago, in a galaxy far far away.
{
	norm.list<-list()
	singles<-0L
	totals<-0L
	others<-0L
    chromosomes<-scanBamHeader(bam.file)[[1]][[1]]
    if (!is.null(restrict)) { chromosomes<-chromosomes[names(chromosomes) %in% restrict] }
	discard <- .processDiscard(discard)

	loose.names <- list()
	for (i in 1:length(chromosomes)) {
		# Filtering out reads.
		chr <- names(chromosomes)[i]
		where<-GRanges(chr, IRanges(1, chromosomes[i]))
		reads<-scanBam(bam.file, param=ScanBamParam(what=c("qname", "flag", "pos", "qwidth", "mapq"), which=where, 
			flag=scanBamFlag(isUnmappedQuery=FALSE,	isDuplicate=ifelse(dedup, FALSE, NA))))[[1]]
		reads <- .discardReads(reads, discard[[chr]])
	    keep<-reads$mapq >= minq & !is.na(reads$mapq) 
		reads$mapq <- NULL
		for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }

		# Processing higher-order statistics.
		totals <- totals + length(reads$flag)
		is.single <- bitwAnd(reads$flag, 0x1)==0L
		singles <- singles + sum(is.single)
 	    for (x in names(reads)) { reads[[x]] <- reads[[x]][!is.single] }	
		okay <- .yieldInterestingBits(reads, where, diag=TRUE)
		norm.list[[i]] <- okay$size
	
		# Storing loose names.
		leftovers <- reads$qname[!okay$is.ok]
		freqs <- table(leftovers)
		others <- others + as.integer(sum(freqs==2L))
		loose.names[[i]] <- names(freqs)[freqs==1L]
	}

	# Cleaning up the loose names; figuring out if they're unmapped or interchromosomal.
	# If there's only 1, then the other read is unmapped somewhere so it's the former. If there's
	# 2, then they must be on different chromosomes to slip through the net.
	loose.names <- unlist(loose.names)
	occurrences <- table(loose.names)
	one.unmapped <- sum(occurrences==1L)
	inter.chr <- sum(occurrences==2L)

	# Returning sizes and some diagnostic data.
    return(list(sizes=unlist(norm.list), diagnostics=c(total=totals, single=singles, unoriented=others,
			mate.unmapped=one.unmapped, inter.chr=inter.chr)))
}

##################################

.extractPET <- function(bam.file, where, dedup=FALSE, minq=0, na.rm=TRUE, discard=NULL)
# A function to extract PET data for a particular chromosome. Synchronisation is expected.
# We avoid sorting by name  as it'd mean we have to process the entire genome at once 
# (can't go chromosome-by-chromosome).  This probably results in increased memory usage 
# across the board, and it doesn't fit in well with the rest of the pipelines which assume 
# coordinate sorting.
# 
# written by Aaron Lun
# 8 December 2013
{
    reads <-scanBam(bam.file, param=ScanBamParam(what=c("qname", "flag", "pos", "qwidth", "mapq"),
		which=where, flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=ifelse(dedup, FALSE, NA), 
		isPaired=TRUE, hasUnmappedMate=FALSE)))[[1]]
    keep<-reads$mapq >= minq
	if (na.rm) { keep<-keep & !is.na(reads$mapq) }
	reads$mapq<-NULL
	for (x in names(reads)) { reads[[x]] <- reads[[x]][keep] }
	reads <- .discardReads(reads, discard)
	.yieldInterestingBits(reads, where)
}

.yieldInterestingBits <- function(reads, where, diag=FALSE)
# This figures out what the interesting reads are, given a list of 
# read names, flags, positions and qwidths. In particular, reads have
# to be properly paired in an inward conformation, without one read
# overrunning the other (if the read lengths are variable).
#
# written by Aaron Lun
# 15 April 2014
{ 
	# Figuring out the strandedness, and matching the reads.
 	should.be.left <- bitwAnd(reads$flag, 0x10) == 0L
	should.be.right <- !should.be.left
	corresponding <- match(reads$qname[should.be.left], reads$qname[should.be.right])
	hasmatch <- !is.na(corresponding)

	# Selecting. Anything that survives the above should be forward in 'left', and reverse in 'right'.
	fpos <- reads$pos[should.be.left][hasmatch]
	fwidth <- reads$qwidth[should.be.left][hasmatch]
	rpos <- reads$pos[should.be.right][corresponding[hasmatch]]
	rwidth <- reads$qwidth[should.be.right][corresponding[hasmatch]]

	# Allowing only valid PETs.
	fend <- pmin(fpos+fwidth, end(where)+1L)
	rend <- pmin(rpos+rwidth, end(where)+1L)
    valid <- fpos <= rpos & fend <= rend
	total.size <- rend - fpos
	fpos <- fpos[valid]
	total.size <- total.size[valid]

	if (diag) { 
		is.ok <- logical(length(reads$flag))
		is.ok[should.be.left][hasmatch][valid] <- TRUE
		is.ok[should.be.right][corresponding[hasmatch]][valid] <- TRUE
		return(list(pos=fpos, size=total.size, is.ok=is.ok))
	}
	return(list(pos=fpos, size=total.size))
}

##################################

.extractBrokenPET <- function(bam.file, where, dedup=FALSE, minq=0, na.rm=TRUE, discard=NULL, use.first=TRUE) 
# A function to extract PET data, but as single-end data (i.e. only using one of the reads).
# Useful when paired-end data has gone completely off the rails.
{
	.extractSET(bam.file, where, dedup=dedup, minq=minq, na.rm=na.rm, discard=discard,
		isPaired=TRUE, isFirstMateRead=use.first, isSecondMateRead=!use.first)
}

##################################
