profileSites <- function(bam.files, regions, range=5000, ext=100, weight=1, 
	param=readParam(), use.strand=TRUE, match.strand=FALSE)
# This is a function to compute the profile around putative binding sites. The 5' edge of the
# binding site is identified by counting reads into a window of size `width`, on the left and
# right of a given position, and determining if the right/left ratio is greater than 5. It then
# records the coverage of the resulting bases, up to `range`.
#
# written by Aaron Lun
# created 2 July 2012
# last modified 10 February 2015
{
	weight <- as.double(weight)
	if(length(weight) != length(regions)) { weight <- rep(weight, length.out=length(regions)) }
	nbam <- length(bam.files)
	paramlist <- .makeParamList(nbam, param)

	# A bit of work for strand-specificity.
	if (match.strand) { 
		use.strand <- TRUE
	    for (i in 1:nbam) { 
			if (length(paramlist[[i]]$forward)) { stop("set forward=NULL in param for strand-specific profiling") } 
		}
	}
	if (use.strand) { 
		reverse <- strand(regions)=="-"
		if (any(reverse)) {
			reverse <- as.logical(reverse)
			rregs <- regions[reverse]
			start(rregs) <- end(rregs) # Using the 5' end of the reverse-stranded region.
			rprof <- Recall(bam.files=bam.files, regions=rregs, range=range, ext=ext, weight=weight[reverse], 
				param=reformList(paramlist, forward=ifelse(match.strand, FALSE, NA)), 
				use.strand=FALSE, match.strand=FALSE)
			if (any(!reverse)) { 
				fprof <- Recall(bam.files=bam.files, regions=regions[!reverse], range=range, ext=ext, weight=weight[!reverse],
 					param=reformList(paramlist, forward=ifelse(match.strand, TRUE, NA)), 
					use.strand=FALSE, match.strand=FALSE)
			} else { 
				fprof <- 0 
			}
			prop.rev <- sum(reverse)/length(reverse)
			return(fprof * (1-prop.rev) + rev(rprof) * prop.rev) # Flipping the profile.
		}
	}

	# Setting up.
	extracted.chrs <- .activeChrs(bam.files, paramlist[[1]]$restrict)
	ext.data <- .collateExt(nbam, ext)
	range <- as.integer(range)
	if (range <= 0L) { stop("range should be positive") }
	total.profile <- 0
	indices <- split(1:length(regions), seqnames(regions))
		
	# Running through the chromosomes.
	for (i in 1:length(extracted.chrs)) {
		chr <- names(extracted.chrs)[i]
		chosen <- indices[[chr]]
		if (!length(chosen)) { next }
		outlen <- extracted.chrs[i]
		where <- GRanges(chr, IRanges(1L, outlen))

        # Reading in the reads for the current chromosome for all the BAM files.
		starts <- ends <- list()
		for (b in 1:nbam) {
			curpar <- paramlist[[b]]
            if (curpar$pe!="both") {
				if (curpar$pe=="none") { 
					reads <- .extractSE(bam.files[b], where=where, param=curpar)
				} else {
					reads <- .extractBrokenPE(bam.files[b], where=where, param=curpar)
				}
   				extended <- .extendSE(reads, ext=ext.data$ext[b])
				start.pos <- extended$start
				end.pos <- extended$end
			} else {
                if (!is.na(curpar$rescue.ext)) {
					out <- .rescuePE(bam.files[b], where=where, param=curpar)
				} else {
					out <- .extractPE(bam.files[b], where=where, param=curpar)
				}
				start.pos <- out$pos
				end.pos <- out$pos + out$size - 1L
			}
			checked <- .checkFragments(start.pos, end.pos, final=ext.data$final, chrlen=outlen)
			start.pos <- checked$start
			end.pos <- checked$end

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
