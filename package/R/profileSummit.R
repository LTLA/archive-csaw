profileSummit <- function(bam.files, ext=100, width=10000, res=50, min.depth=20, param=readParam()) 
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

	min.depth <- as.integer(min.depth)
	res <- as.integer(res)
	actual.width <- as.integer(width/res)
	ext <- as.integer(ext)
    if (res<=0) { stop("bin size resolution must be positive") }
    if (actual.width <=0) { stop("smoothing width must be positive") }
    if (ext <=0) { stop("extension length must be positive") }

	blen <- length(bam.files)
	total.profile <- total.freq <- 0

	# Running through the chromosomes.
    for (i in 1:length(extracted$chrs)) {
		chr <- names(extracted$chrs)[i]
		outlen <- extracted$chrs[i]
		where <- GRanges(chr, IRanges(1L, outlen))

		total.pts <- as.integer((outlen-1)/res)+1L
		total.cov <- integer(total.pts)

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
				positions <- as.integer(start.pos + ext/2)
			} else {
                if (rescue.pairs) {
					out <- .rescuePET(bam.files[bf], where=where, dedup=dedup, minq=minq,
						max.frag=max.frag, ext=rescue.ext, discard=extracted$discard[[chr]])
				} else {
					out <- .extractPET(bam.files[bf], where=where, dedup=dedup, minq=minq,
						discard=extracted$discard[[chr]], max.frag=max.frag)
				}
				positions <- as.integer(out$pos + out$size/2)
			}
			if (!length(positions)) { next }

			# Setting up the coverage track on each bin (pooling counts across files).
			out <- .Call(cxx_get_rle_counts, positions-res+1L, positions, total.pts, res, TRUE)
			if (is.character(out)) { stop(out) }
			total.cov <- total.cov + out
		}
		
	    # We call the C++ function to aggregate profiles.
		cur.profile <- .Call(cxx_get_profile, total.cov, actual.width, min.depth)
		if (is.character(cur.profile)) { stop(cur.profile) }
		total.profile <- total.profile + cur.profile[[1]]
		total.freq <- total.freq + cur.profile[[2]]
    }

	print(total.freq)
	print(total.profile)
	# Cleaning up and returning the profiles.
    return(list(span=1:actual.width*res, coverage=(total.profile/total.freq)))
}


