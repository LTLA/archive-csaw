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
                if (.rescueMe(curpar)) { 
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

wwhm <- function(profile, regions, ext=100, param=readParam(), proportion=0.5, len=NULL)
# This function computes the window width at half its maximum. This uses
# the output of profileSites to get the full width of the peak; it then
# subtracts twice the extension length to obtain the window width. 
# 
# written by Aaron Lun
# created 2 March 2015
{
	if (proportion <= 0 | proportion >= 1) { stop("proportion should be between 0 and 1") }
	is.max <- which.max(profile)
	cutoff <- proportion * profile[is.max]
	above.max <- profile >= cutoff

	# Getting the width of the peak at half-max.
	out <- rle(above.max)
	ends <- cumsum(out$lengths)
	starts <- c(1L, ends[-length(ends)]+1L)
	for (ok in which(out$values)) { 
		if (is.max <= ends[ok] & is.max >= starts[ok]) {
			chosen.start <- starts[ok]
			chosen.end <- ends[ok]
			break
		}
	}
	if (chosen.end==length(profile) || chosen.start==1L) { 
		warning("width at specified proportion exceeds length of profile")
	}
	peak.width <- chosen.end - chosen.start + 1L

	# Getting the median size of the regions.
	ref.size <- median(width(regions))

	# To get the average extension length across libraries, via getWidths.
	nbam <- max(length(ext), ifelse(is.list(param), length(param), 1))
	paramlist <- .makeParamList(nbam, param)
	dim(paramlist) <- c(nbam, 1)
	colnames(paramlist) <- "param"
	ext.data <- .collateExt(nbam, ext)
	dummy.data <- SummarizedExperiment(matrix(0, ncol=length(ext), nrow=1),
		rowData=GRanges("chrA", IRanges(1, 1)), 
		colData=DataFrame(ext=ext.data$ext, paramlist),
		exptData=List(final.ext=ext.data$final))
	ext.len <- getWidths(dummy.data, len=len)

	# Computing the window size. Add 2 to ensure one base overlaps the 
	# most extreme fragments on both sides. Subtract ref.size-1 as random
	# distribution of summits within maximal windows widens the peak.
	max(peak.width - ext.len*2L + 2L - ref.size + 1L, 1L)
}
