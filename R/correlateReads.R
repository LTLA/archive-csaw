correlateReads <- function(bam.files, max.dist=1000, cross=TRUE, param=readParam()) 
# This is just a function to calculate the autocorrelation between reads of different strands (or
# between reads in general). Note that the BAM files must be sorted. It will calculate the values 
# required for computation of the correlation function across all chromosomes, then it will crunch all 
# the data together to get the correlation coefficient at each distance. I haven't used R's native 
# acf/ccf() function because it doesn't handle the large inputs from BAM efficiently.
#
# written by Aaron Lun
# created 2 July 2012
# last modified 5 August 2016
{
	nbam <- length(bam.files)
    if (is.list(param)) { 
        .Deprecated(msg="supplying a list of readParam objects is deprecated, using first element only")
        param <- param[[1]]
    }
	extracted.chrs <- .activeChrs(bam.files, param$restrict)

	max.dist <- as.integer(max.dist)
	if (max.dist <=0) { stop("maximum distance must be positive") }
	total.cor <- numeric(max.dist+1L)
	total.read.num <- 0L

	for (i in seq_along(extracted.chrs)) {
		chr <- names(extracted.chrs)[i]
		where <- GRanges(chr, IRanges(1L, extracted.chrs[i]))
        total.len <- extracted.chrs[i] + 1L # Length of the conceptual vector to compute the correlations.
        if (total.len < 2L) { next } # No way to compute variance if the vector's too small.

		# Reading in the reads for the current chromosome for all the BAM files.
        bp.out <- bplapply(seq_len(nbam), FUN=.correlate_reads,
                           bam.files=bam.files, where=where, param=param, 
                           total.len=total.len,
                           BPPARAM=param$BPPARAM)

        all.f <- lapply(bp.out, "[[", "forward")
        all.r <- lapply(bp.out, "[[", "reverse")
        if (!cross) {
            all.f <- mapply("c", all.f, all.r, SIMPLIFY=FALSE)
        }
        forward.reads <- sum(sapply(bp.out, "[[", "nforward"))
        num.reads <- forward.reads + sum(sapply(bp.out, "[[", "nreverse"))

		# Assembling RLEs (with some protection from empties). We need reads for any correlation and 
		# reads on both strands to get cross-correlations. If we're doing cross-correlations, then
		# we compare between strands; if we're doing autocorrelations, we compare within all reads.		
		if (num.reads==0L) { next; }
		if (cross && (forward.reads==0L || forward.reads==num.reads)) { next } # correlations undefined, so they don't contribute to total.read.num.
		all.f <- rle(sort(do.call(c, all.f)))
		if (cross) {
			all.r <- rle(sort(do.call(c, all.r)))
		} else {
			all.r <- all.f
		}

		# We call the C++ function to compute correlations. 
		ccfs <- .Call(cxx_correlate_reads, all.f$values, all.f$lengths, all.r$values, all.r$lengths, max.dist, total.len)
		if (is.character(ccfs)) { stop(ccfs) }

		# Returning some output. Note that the coefficient is weighted according to the number
		# of reads on each chromosome, as described in as described by Kharchenko et al. (2008).
		total.read.num <- total.read.num+num.reads
		total.cor <- total.cor+ccfs*num.reads
	}

	# Cleaning up and returning the correlations.
	if (total.read.num) { total.cor <- total.cor/total.read.num }
	return(total.cor)
}

.correlate_reads <- function(bf, bam.files, where, param, total.len) {
    if (param$pe=="both") {
        reads <- .getPairedEnd(bam.files[bf], where=where, param=param, with.reads=TRUE)
    } else {
        reads <- .getSingleEnd(bam.files[bf], where=where, param=param)
    }
    
    forward.pos <- reads$forward$pos
    forward.pos[forward.pos < 1L] <- 1L
    reverse.pos <- reads$reverse$pos + reads$reverse$qwidth
    reverse.pos[reverse.pos > total.len] <- total.len

    num.reads <- length(forward.pos)+length(reverse.pos)
    forward.reads <- length(forward.pos)
    return(list(forward=forward.pos, reverse=reverse.pos, 
                nforward=length(forward.pos), nreverse=length(reverse.pos)))
}

