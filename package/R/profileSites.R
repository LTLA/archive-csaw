profileSites <- function(bam.files, regions, range=5000, ext=100, weight=1, 
    final.ext=NULL, param=readParam()) 
# This is a function to compute the profile around putative binding sites. The 5' edge of the
# binding site is identified by counting reads into a window of size `width`, on the left and
# right of a given position, and determining if the right/left ratio is greater than 5. It then
# records the coverage of the resulting bases, up to `max.dist`.
#
# written by Aaron Lun
# created 2 July 2012
# last modified 13 December 2014
{
	nbam <- length(bam.files)
	paramlist <- .makeParamList(nbam, param)
	extracted.chrs <- .activeChrs(bam.files, paramlist[[1]]$restrict)
	ext.data <- .collateExt(nbam, ext, final.ext)

	# Splitting up the regions.
	indices <- split(1:length(regions), seqnames(regions))
	weight <- as.double(weight)
	if(length(weight) != length(regions)) { weight <- rep(weight, length.out=length(regions)) }
	
	range <- as.integer(range)
	if (range <= 0L) { stop("range should be positive") }
	blen <- length(bam.files)
	total.profile <- 0
		
	# Running through the chromosomes.
    for (i in 1:length(extracted.chrs)) {
		chr <- names(extracted.chrs)[i]
		chosen <- indices[[chr]]
		if (!length(chosen)) { next }
		outlen <- extracted.chrs[i]
		where <- GRanges(chr, IRanges(1L, outlen))

        # Reading in the reads for the current chromosome for all the BAM files.
		starts <- ends <- list()
		for (b in 1:blen) {
			curpar <- paramlist[[b]]
            if (curpar$pe!="both") {
				if (curpar$pe=="none") { 
					reads <- .extractSE(bam.files[b], where=where, param=curpar)
				} else {
					reads <- .extractBrokenPE(bam.files[b], where=where, param=curpar)
				}
   				extended <- .extendSE(reads, chrlen=outlen, ext.info=ext.data[b,])
				start.pos <- extended$start
				end.pos <- extended$end
			} else {
                if (curpar$rescue.pairs) {
					out <- .rescuePE(bam.files[b], where=where, param=curpar)
				} else {
					out <- .extractPE(bam.files[b], where=where, param=curpar)
				}
				start.pos <- out$pos
				end.pos <- out$pos + out$size - 1L
			}

			if (!length(start.pos)) { next }
			ix <- length(starts) + 1L
			starts[[ix]] <- pmax(start.pos, 1L) # Avoid considering off ends of chromosomes.
 			ends[[ix]] <- pmin(end.pos, outlen)
		}
			
		# Pulling out the regions.
		all.starts <- start(regions)[chosen]
		os <- order(all.starts)
		all.starts <- all.starts[os]
		all.weights <- weight[chosen][os]

	    # We call the C++ functions to aggregate profiles.
		starts <- unlist(starts)
		ends <- unlist(ends)
		if (!length(starts)) { next }
		cur.profile <- .Call(cxx_get_profile, starts, ends, all.starts, all.weights, range) 
		if (is.character(cur.profile)) { stop(cur.profile) }
		total.profile <- total.profile + cur.profile
    }

	# Cleaning up and returning the profiles. We divide by 2 to get the coverage,
	# as total.profile counts both sides of each summit (and is twice as large as it should be).
	out <- total.profile/length(regions)
	names(out) <- (-range):range
    return(out)
}
