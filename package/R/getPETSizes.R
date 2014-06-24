getPETSizes <- function(bam.file, dedup=FALSE, minq=NULL, restrict=NULL, discard=NULL) 
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
	stopifnot(length(bam.file)==1L)
    extracted <- .processIncoming(bam.file, restrict, discard)
	loose.names.1 <- list()
	loose.names.2 <- list()
	one.unmapped <- 0L

	for (i in 1:length(extracted$chrs)) {
		chr <- names(extracted$chrs)[i]
		where <- GRanges(chr, IRanges(1L, extracted$chrs[i]))
		reads <- .extractSET(bam.file, extras=c("qname", "flag", "isize"), where=where, 
			dedup=dedup, discard=extracted$discard[[chr]], minq=minq)

		# Getting rid of unpaired reads.
		totals <- totals + length(reads$flag)
		is.single <- bitwAnd(reads$flag, 0x1)==0L
		singles <- singles + sum(is.single)
 		only.one <- bitwAnd(reads$flag, 0x8)!=0L
		one.unmapped <- one.unmapped + sum(!is.single & only.one)
	    for (x in names(reads)) { reads[[x]] <- reads[[x]][!is.single & !only.one] }	
		
		# Identifying valid reads.
		okay <- .yieldInterestingBits(reads, extracted$chrs[i], diag=TRUE)
		norm.list[[i]] <- okay$size
		left.names <- reads$qname[!okay$is.ok]
		left.flags <- reads$flag[!okay$is.ok]
		
		# Setting up some more filters.
		on.same.chr <- reads$isize[!okay$is.ok]!=0L
		is.first <- bitwAnd(left.flags, 0x40)!=0L
		is.second <- bitwAnd(left.flags, 0x80)!=0L
		
		# Identifying improperly orientated pairs and reads with unmapped counterparts.
		leftovers.first <- left.names[on.same.chr & is.first]
		leftovers.second <- left.names[on.same.chr & is.second]
		has.pair <- sum(leftovers.first %in% leftovers.second)
		others <- others + has.pair
		one.unmapped <- one.unmapped + length(leftovers.first) + length(leftovers.second) - 2L*has.pair 

		# Collecting the rest to match inter-chromosomals.
		loose.names.1[[i]] <- left.names[!on.same.chr & is.first]
		loose.names.2[[i]] <- left.names[!on.same.chr & is.second]
	}

	# Checking whether a read is positively matched to a mapped counterpart on another chromosome.
	# If not, then it's just an read in an unmapped pair.
	loose.names.1 <- unlist(loose.names.1)
	loose.names.2 <- unlist(loose.names.2)
	inter.chr <- sum(loose.names.1 %in% loose.names.2)
	one.unmapped <- one.unmapped + length(loose.names.2) + length(loose.names.1) - inter.chr*2L

	# Returning sizes and some diagnostic data.
    return(list(sizes=unlist(norm.list), diagnostics=c(total=totals, single=singles, unoriented=others,
			mate.unmapped=one.unmapped, inter.chr=inter.chr)))
}

##################################

.extractPET <- function(bam.file, where, dedup, minq, max.frag=Inf, discard=NULL)
# A function to extract PET data for a particular chromosome. Synchronisation is expected.
# We avoid sorting by name  as it'd mean we have to process the entire genome at once 
# (can't go chromosome-by-chromosome).  This probably results in increased memory usage 
# across the board, and it doesn't fit in well with the rest of the pipelines which assume 
# coordinate sorting.
# 
# written by Aaron Lun
# 8 December 2013
{
	reads <- .extractSET(bam.file, extras=c("qname", "flag"), where=where, dedup=dedup,
		minq=minq, isPaired=TRUE, hasUnmappedMate=FALSE, discard=discard)
	.yieldInterestingBits(reads, max(end(where)), max.frag=max.frag)
}

.yieldInterestingBits <- function(reads, clen, max.frag=Inf, diag=FALSE)
# This figures out what the interesting reads are, given a list of 
# read names, flags, positions and qwidths. In particular, reads have
# to be properly paired in an inward conformation, without one read
# overrunning the other (if the read lengths are variable).
#
# written by Aaron Lun
# 15 April 2014
{ 
	if (max.frag <= 0L) {
		return(list(pos=integer(0), size=integer(0), is.ok=logical(length(reads$flag))))
	}

	# Figuring out the strandedness and the read number.
 	is.forward <- bitwAnd(reads$flag, 0x10) == 0L 
 	is.mate.reverse <- bitwAnd(reads$flag, 0x20) != 0L
	should.be.left <- is.forward & is.mate.reverse
	should.be.right <- !is.forward & !is.mate.reverse
    is.first <- bitwAnd(reads$flag, 0x40) != 0L
	is.second <- bitwAnd(reads$flag, 0x80) != 0L	
	stopifnot(all(is.first!=is.second))

	# Matching the reads in each pair so only valid PETs are formed.
	set.first.A <- should.be.left & is.first
	set.second.A <- should.be.right & is.second
	set.first.B <- should.be.left & is.second
	set.second.B <- should.be.right & is.first
	corresponding.1 <- match(reads$qname[set.first.A], reads$qname[set.second.A])
	corresponding.2 <- match(reads$qname[set.first.B], reads$qname[set.second.B])

	hasmatch <- logical(length(reads$flag))
	hasmatch[set.first.A] <- !is.na(corresponding.1)
	hasmatch[set.first.B] <- !is.na(corresponding.2)
	corresponding <- integer(length(reads$flag))
	corresponding[set.first.A] <- which(set.second.A)[corresponding.1]
	corresponding[set.first.B] <- which(set.second.B)[corresponding.2]

	# Selecting the read pairs.
	fpos <- reads$pos[hasmatch]
	fwidth <- reads$qwidth[hasmatch]
	rpos <- reads$pos[corresponding[hasmatch]]
	rwidth <- reads$qwidth[corresponding[hasmatch]]

	# Allowing only valid PETs.
	fend <- pmin(fpos+fwidth, clen+1L)
	rend <- pmin(rpos+rwidth, clen+1L)
	total.size <- rend - fpos
    valid <- fpos <= rpos & fend <= rend & total.size <= max.frag
	fpos <- fpos[valid]
	total.size <- total.size[valid]

	if (diag) { 
		is.ok <- logical(length(reads$flag))
		is.ok[hasmatch][valid] <- TRUE
		is.ok[corresponding[hasmatch]][valid] <- TRUE
		return(list(pos=fpos, size=total.size, is.ok=is.ok))
	}
	return(list(pos=fpos, size=total.size))
}

##################################

.extractBrokenPET <- function(bam.file, where, dedup, minq, discard=NULL, use.first=TRUE) 
# A function to extract PET data, but as single-end data (i.e. only using one
# of the reads).  Useful when paired-end data has gone completely off the
# rails.
{
	.extractSET(bam.file, where=where, dedup=dedup, minq=minq, discard=discard,
		isPaired=TRUE, isFirstMateRead=use.first, isSecondMateRead=!use.first)
}

.rescuePET <- function(bam.file, where, dedup, minq, discard=NULL, ext=100, max.frag=Inf) 
# A function to extract PET data where possible, but to rescue those that
# are invalid by using them as single-end data with read extension. Those
# reads that form invalid pairs are broken up and the read with the better
# MAPQ is chosen. Any single read (due to filtering or whatever) is used as-is.
# Interchromosomal pairs get counted once on each chromosome.
#
# written by Aaron Lun
# 13 May, 2014
{
	reads <- .extractSET(bam.file, extras=c("qname", "flag", "mapq"), where=where, dedup=dedup,
		minq=minq, isPaired=TRUE, discard=discard)
	output <- .yieldInterestingBits(reads, max(end(where)), diag=TRUE, max.frag=max.frag)

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

	# Returning the loot. Yar!	
	return( list( pos=c(output$pos, ifelse(bitwAnd(reads$flag[additor], 0x10)==0L, reads$pos[additor], 
					reads$pos[additor]+reads$qwidth[additor]-ext)),
		      size=c(output$size, rep(ext, sum(additor))) ) )
}

##################################
