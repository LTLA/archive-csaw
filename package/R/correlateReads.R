correlateReads <- function(bam.files, max.dist=1000, dedup=FALSE, minq=0, cross=TRUE, restrict=NULL, 
	discard=NULL, pet=c("none", "first", "second")) 
# This is just a function to calculate the autocorrelation between reads of different strands (or
# between reads in general). Note that the BAM files must be sorted. It will calculate the values 
# required for computation of the correlation function across all chromosomes, then it will crunch all 
# the data together to get the correlation coefficient at each distance. I haven't used R's native 
# acf/ccf() function because it doesn't handle the large inputs from BAM efficiently.
#
# written by Aaron Lun
{
    extracted <- .processIncoming(bam.files, restrict, discard)
	max.dist<-as.integer(max.dist)
    if (max.dist <=0 ) { stop("maximum distance must be positive") }
    total.cor <- numeric(max.dist+1L)
    total.read.num <- 0L;
	pet <- match.arg(pet)

    for (i in 1:length(extracted$chrs)) {
		if (extracted$chrs[i]<2L) { next } # No way to compute variance if there's only one base.
		chr <- names(extracted$chrs)[i]
		where <- GRanges(chr, IRanges(1L, extracted$chrs[i]))

        # Reading in the reads for the current chromosome for all the BAM files.
		all.f<-all.r<-list()
		num.reads<-0
		forward.reads<-0
		for (b in 1:length(bam.files)) {
			if (pet=="none") { 
				reads<-.extractSET(bam.files[b], where=where, dedup=dedup, minq=minq, discard=extracted$discard[[chr]])
			} else {
				reads<-.extractBrokenPET(bam.files[b], where=where, dedup=dedup, minq=minq, discard=extracted$discard[[chr]],
					use.first=(pet=="first"))
			}
			forwards<-reads$strand=="+"
			num.reads<-num.reads+length(forwards)
			forward.reads<-forward.reads+sum(forwards)
			if (!all(forwards)) { reads$pos[!forwards]<-pmin(extracted$chrs[i], reads$pos[!forwards]+reads$qwidth[!forwards]) }
			if (cross) {
				all.f[[b]]<-reads$pos[forwards]
				all.r[[b]]<-reads$pos[!forwards]
			} else { all.f[[b]]<-reads$pos }
		}

		# Assembling RLEs (with some protection from empties). We need reads for any correlation and 
		# reads on both strands to get cross-correlations. If we're doing cross-correlations, then
		# we compare between strands; if we're doing autocorrelations, we compare within all reads.		
		if (num.reads==0) { next; }
		if (cross && (forward.reads==0 || forward.reads==num.reads)) { next }
		all.f<-rle(sort(do.call(c, all.f)))
		if (cross) {
			all.r<-rle(sort(do.call(c, all.r)))
		} else {
			all.r<-all.f
		}

	    # We call the C++ function to compute correlations. 
		ccfs<-.Call(cxx_correlate_reads, all.f$values, all.f$lengths, all.r$values, all.r$lengths, max.dist, extracted$chrs[i])
		if (is.character(ccfs)) { stop(ccfs) }

		# Returning some output. Note that the coefficient is weighted according to the number
        # of reads on each chromosome, as described in as described by Kharchenko et al. (2008).
        total.read.num<-total.read.num+num.reads;
        total.cor<-total.cor+ccfs*num.reads;
    }

	# Cleaning up and returning the correlations.
	if (total.read.num) { total.cor<-total.cor/total.read.num }
    return(total.cor)
}


