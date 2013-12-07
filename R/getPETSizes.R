getPETSizes <- function(bam.file, dedup=FALSE, minq=0, restrict=NULL) 
# This function takes a BAM file and reads it to parse the size of the PET fragments. It then
# returns a vector of sizes which can be plotted for diagnostics. The length of the vector
# will also tell you how many read pairs were considered valid. The total number of reads, the
# number of singletons and the number of interchromosomal pairs is also reported.
# 
# written by Aaron Lun
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
	
		# Pulling out other diagnostics.
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

.extractPET <- function(bam.file, where, dedup=FALSE, minq=0) 
# A function to extract PET data for a particular chromosome. We rely on 
# properly synchronised read pair data for convenience. Otherwise, we'd
# have to figure out which reads were paired. This becomes inconveniet
# if sorted by position; and, if sorted by name, it'd mean we have to process
# the entire genome at once (can't go chromosoem-by-chromosome).
{
    reads<-scanBam(bam.file, param=ScanBamParam(what=c("pos", "mpos", "isize", "mapq"),
		which=where, flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=ifelse(dedup, FALSE, NA), 
		isPaired=TRUE, hasUnmappedMate=FALSE, isMinusStrand=FALSE, isMateMinusStrand=TRUE)))[[1]]

	# Protecting against undefined behavior when fragment=read length. Which one is the 
	# rightmost or leftmost segment may not be defined as they only differ in strandedness.
	# Unfortunately, it won't filter out scenarios where the reverse read is nested in the forward read.
	# This should be only possible when working with variable length sequencing data (e.g. Roche).
	reads$isize<-abs(reads$isize)

	invalid<-reads$pos > reads$mpos 	# removes reverse reads which lie before the forward read
	inter.chr<-reads$isize==0L			# gets rid of reads from different chromosomes
	good.q<-reads$mapq >= minq 			# gets rid of pairs with low mapping quality (assuming they are equal between reads in a pair) 
	keep<-!(invalid | inter.chr) & good.q
	return(list(pos=reads$pos[keep], size=reads$isize[keep]))
}
