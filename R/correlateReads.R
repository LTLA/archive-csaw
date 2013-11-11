correlateReads <- function(bam.files, max.dist=1000, dedup=FALSE, minq=0, cross=TRUE, restrict=NULL) 
# This is just a function to calculate the autocorrelation between reads of different strands (or
# between reads in general). Note that the BAM files must be sorted. It will calculate the values 
# required for computation of the correlation function across all chromosomes, then it will crunch all 
# the data together to get the correlation coefficient at each distance. I haven't used R's native 
# acf/ccf() function because it doesn't handle the large inputs from BAM efficiently.
#
# written by Aaron Lun
{
	# Loading various details of the biological system.
	bf<-lapply(bam.files, FUN=function(b) { open(BamFile(b)); })
    chromosomes<-scanBamHeader(bam.files[1])[[1]][[1]]
	if (!is.null(restrict)) { chromosomes<-chromosomes[names(chromosomes) %in% restrict] }
	
	# Loading containers.
	max.dist<-as.integer(max.dist)
    if (max.dist <=0 ) { stop("maximum distance must be positive"); }
    total_cor<-numeric(max.dist+1L);
    total_read_num<-0L;

    for (i in 1:length(chromosomes)) {
		if (chromosomes[i]==0) { next }

        # Reading in the reads for the current chromosome for all the BAM files.
		my.chr<-GRanges(names(chromosomes)[i], IRanges(1, chromosomes[i]))
		all.f<-all.r<-list()
		num.reads<-0
		forward.reads<-0
		for (b in 1:length(bf)) {
			reads<-.extractSET(bf[[b]], where=my.chr, dedup=dedup, minq=minq)
			forwards<-(reads$strand=="+");
			num.reads<-num.reads+length(forwards)
			forward.reads<-forward.reads+sum(forwards)
			if (!all(forwards)) { reads$pos[!forwards]<-pmin(chromosomes[i], reads$pos[!forwards]+reads$qwidth[!forwards]);	}
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
		ccfs<-.Call("R_correlate_reads", all.f$values, all.f$lengths, all.r$values, all.r$lengths, max.dist, 
			chromosomes[i], PACKAGE="csaw");
		if (is.character(ccfs)) { stop(ccfs); }

		# Returning some output. Note that the coefficient is weighted according to the number
        # of reads on each chromosome, as described in as described by Kharchenko et al. (2008).
        total_read_num<-total_read_num+num.reads;
        total_cor<-total_cor+ccfs*num.reads;
    }

	# Cleaning up and returning the correlations.
	for (b in bf) { close(b) }
	if (total_read_num) { total_cor<-total_cor/total_read_num }
    return(total_cor)
}

###################################################################################################
# End.
