profileSummit <- function(bam.files, ext=100, width=5000, res=50, min.depth=1, param=readParam()) 
# This is a function to compute the profile around putative binding sites. The 5' edge of the
# binding site is identified by counting reads into a window of size `width`, on the left and
# right of a given position, and determining if the right/left ratio is greater than 5. It then
# records the coverage of the resulting bases, up to `max.dist`.
#
# written by Aaron Lun
# 2 July 2012
# modified 1 September, 2014
{
    extracted <- .processIncoming(bam.files, param$restrict, param$discard)
	pet <- param$pet
	minq <- param$minq
	dedup <- param$dedup
	rescue.pairs <- param$rescue.pairs
	rescue.ext <- param$rescue.ext
	max.frag <- param$max.frag

	res <- as.integer(res)
	actual.width <- as.integer(width/res)
	ext <- as.integer(ext)
    if (res<=0) { stop("bin size resolution must be positive") }
    if (actual.width <=0) { stop("smoothing width must be positive") }
    if (ext <=0) { stop("extension length must be positive") }

	min.depth <- as.integer(ceiling(min.depth*res))
	blen <- length(bam.files)
	total.profile <- total.freq <- 0

	# Running through the chromosomes.
    for (i in 1:length(extracted$chrs)) {
		chr <- names(extracted$chrs)[i]
		outlen <- extracted$chrs[i]
		where <- GRanges(chr, IRanges(1L, outlen))

		total.pts <- as.integer((outlen-1)/res)+1L
		total.cov <- integer(total.pts)
		starts <- ends <- list()

        # Reading in the reads for the current chromosome for all the BAM files.
		for (b in 1:blen) {
            if (pet!="both") {
				if (pet=="none") { 
					reads <- .extractSET(bam.files[b], where=where, dedup=dedup, minq=minq, 
						discard=extracted$discard[[chr]])
				} else {
					reads <- .extractBrokenPET(bam.files[b], where=where, dedup=dedup, minq=minq, 
						discard=extracted$discard[[chr]], use.first=(pet=="first"))
				} 
				start.pos <- ifelse(reads$strand=="+", reads$pos, reads$pos + reads$qwidth - ext)
				end.pos <- start.pos + ext - 1L
			} else {
                if (rescue.pairs) {
					out <- .rescuePET(bam.files[b], where=where, dedup=dedup, minq=minq,
						max.frag=max.frag, ext=rescue.ext, discard=extracted$discard[[chr]])
				} else {
					out <- .extractPET(bam.files[b], where=where, dedup=dedup, minq=minq,
						discard=extracted$discard[[chr]], max.frag=max.frag)
				}
				start.pos <- out$pos
				end.pos <- out$pos + out$size - 1L
			}

			if (!length(start.pos)) { next }
			ix <- length(starts) + 1L
			start.pos <- start.pos - res + 1L
			starts[[ix]] <- start.pos 
 			ends[[ix]] <- end.pos
		}
		
	    # We call the C++ functions to aggregate profiles.
		starts <- unlist(starts)
		ends <- unlist(ends)
		if (!length(starts)) { next }
		cur.profile <- .Call(cxx_get_profile, starts, ends, total.pts, res, actual.width, min.depth)
		if (is.character(cur.profile)) { stop(cur.profile) }
		total.profile <- total.profile + cur.profile[[1]]
		total.freq <- total.freq + cur.profile[[2]]
    }

	# Cleaning up and returning the profiles. We multiply by 2 to get the span, as 
	# we need to consider both sides of the summit. We multiply by 2 to get the coverage,
	# as total.freq counts both sides of each summit (and is twice as large as it should be).
    return(list(span=1:actual.width*res*2 + res, coverage=total.profile/total.freq*2))
}


