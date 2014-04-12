getPETSizes <- function(bam.file, dedup=FALSE, minq=0, restrict=NULL) 
# This function takes a BAM file and reads it to parse the size of the PET fragments. It then
# returns a vector of sizes which can be plotted for diagnostics. The length of the vector
# will also tell you how many read pairs were considered valid. The total number of reads, the
# number of singletons and the number of interchromosomal pairs is also reported.
# 
# written by Aaron Lun
# a long long time ago, in a galaxy far far away.
{
	norm.list<-list()
	inter.chr<-0L
	singles<-0L
	one.unmapped<-0L
	totals<-0L
    chromosomes<-scanBamHeader(bam.file)[[1]][[1]]
    if (!is.null(restrict)) { chromosomes<-chromosomes[names(chromosomes) %in% restrict] }

	for (i in 1:length(chromosomes)) {
		where<-GRanges(names(chromosomes)[i], IRanges(1, chromosomes[i]))
		out<-.extractPET(bam.file, where, dedup=dedup, minq=minq)
		norm.list[[i]]<-out$size
	
		# Pulling out other diagnostics, including less-good reads which were left out in .extractPET.
		reads<-scanBam(bam.file, param=ScanBamParam(what=c("flag", "isize", "mapq"), which=where, 
			flag=scanBamFlag(isUnmappedQuery=FALSE,	isDuplicate=ifelse(dedup, FALSE, NA))))[[1]]
		curf<-reads$flag
		paired<-bitwAnd(curf, 0x1)!=0
		curf<-curf[paired]
		both.mapped<-bitwAnd(curf, 0x8)==0
		curf<-curf[both.mapped]
 	   	is.first<-bitwAnd(curf, 0x40)!=0

		# Collating them as necessary.
		totals<-totals+length(reads$flag)
		singles<-singles+sum(!paired)
		one.unmapped<-one.unmapped + sum(!both.mapped)				
		inter.chr<-inter.chr+sum(reads$isize[paired][both.mapped][is.first]==0L)
	}

	# Returning sizes and some diagnostic data.
    return(list(sizes=unlist(norm.list), diagnostics=c(total=totals, single=singles,
			mate.unmapped=one.unmapped, inter.chr=inter.chr)));
}

.extractPET <- function(bam.file, where, dedup=FALSE, minq=0, na.rm=TRUE)
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
	reads$qname <- reads$qname[keep]
	reads$flag <- reads$flag[keep]
	reads$pos<-reads$pos[keep]
	reads$qwidth<-reads$qwidth[keep]
	reads$mapq<-NULL

	# Figuring out the strandedness and read'ness.
 	is.forward <- bitwAnd(reads$flag, 0x10) == 0L
	has.rev.mate <- bitwAnd(reads$flag, 0x20) != 0L
	should.be.left <- is.forward & has.rev.mate
	should.be.right <- !is.forward & !has.rev.mate

	# Matching the two.
	corresponding <- match(reads$qname[should.be.left], reads$qname[should.be.right])
	hasmatch <- !is.na(corresponding)
	fpos <- reads$pos[should.be.left][hasmatch]
	fwidth <- reads$qwidth[should.be.left][hasmatch]
	rpos <- reads$pos[should.be.right][corresponding[hasmatch]]
	rwidth <- reads$qwidth[should.be.right][corresponding[hasmatch]]

	# Allowing only valid PETs.
	fend <- pmin(fpos+fwidth, end(where)+1L)
	rend <- pmin(rpos+rwidth, end(where)+1L)
    valid <- fpos <= rpos & fend <= rend
	total.size <- rend - fpos 
	return(list(pos=fpos[valid], size=total.size[valid]))
}

.extractBrokenPET <- function(bam.file, where, dedup=FALSE, minq=0, na.rm=TRUE, use.first=TRUE) 
# A function to extract PET data, but as single-end data (i.e. only using one of the reads).
# Useful when paired-end data has gone completely off the rails.
{
	.extractSET(bam.file, where, dedup=dedup, minq=minq, na.rm=na.rm, 
		isPaired=TRUE, isFirstMateRead=use.first, isSecondMateRead=!use.first)
}
