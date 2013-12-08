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
		one.unmapped<-one.unmapped<-sum(!both.mapped)				
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
    read1 <-scanBam(bam.file, param=ScanBamParam(what=c("qname", "pos", "qwidth", "mapq"),
		which=where, flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=ifelse(dedup, FALSE, NA), 
		isPaired=TRUE, hasUnmappedMate=FALSE, isMinusStrand=FALSE, isMateMinusStrand=TRUE)))[[1]]
    keep<-read1$mapq >= minq
	if (na.rm) { keep<-keep & !is.na(read1$mapq) }
	read1$qname <- read1$qname[keep]
	read1$pos<-read1$pos[keep]
	read1$qwidth<-read1$qwidth[keep]
	read1$mapq<-NULL

    read2 <-scanBam(bam.file, param=ScanBamParam(what=c("qname", "pos", "qwidth", "mapq"),
		which=where, flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=ifelse(dedup, FALSE, NA), 
		isPaired=TRUE, hasUnmappedMate=FALSE, isMinusStrand=TRUE, isMateMinusStrand=FALSE)))[[1]]
    keep<-read2$mapq >= minq
	if (na.rm) { keep<-keep & !is.na(read2$mapq) }
	read2$qname <- read2$qname[keep]
	read2$strand<-read2$strand[keep]
	read2$pos<-read2$pos[keep]
	read2$qwidth<-read2$qwidth[keep]
	read2$mapq<-NULL

	# Matching the two.
	corresponding <- match(read1$qname, read2$qname)
	hasmatch <- !is.na(corresponding)
	fpos <- read1$pos[hasmatch]
	fwidth <- read1$qwidth[hasmatch]
	rpos <- read2$pos[corresponding[hasmatch]]
	rwidth <- read2$qwidth[corresponding[hasmatch]]

	# Allowing only valid PETs.
	fend <- pmin(fpos+fwidth, end(where)+1L)
	rend <- pmin(rpos+rwidth, end(where)+1L)
    valid <- fpos <= rpos & fend <= rend
	total.size <- rend - fpos 
	return(list(pos=fpos[valid], size=total.size[valid]))
}
