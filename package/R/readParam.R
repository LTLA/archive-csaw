# This defines the readParam class, which specifies the universal
# parameters for read loading in csaw. A class is used here so that
# there are continuous validity checks on the list values.

setClass("readParam", representation(pe="character", 
	max.frag="integer", rescue.pairs="logical", rescue.ext="integer", 
	dedup="logical", minq="integer", forward="logical", 
	restrict="character", discard="GRanges"))

setValidity("readParam", function(object) {
    if (length(object@pe)!=1L || ! object@pe %in%	c("none", "both", "first", "second")) { 
		return("PE specification must be a character scalar of 'none', 'both', 'first' or 'second'") 
	}
   	if (length(object@max.frag)!=1L || object@max.frag <= 0L) {
		return("maximum fragment specifier must be a positive integer")
	} 

	if (length(object@rescue.pairs)!=1L) {
		return("pair rescue specification must be a logical scalar")
	}
	if (length(object@rescue.ext)!=1L) {
		if (is.na(object@rescue.ext)) { 
   			if (object@rescue.pairs) { 
				return("extension length must be specified for rescuing")
			}
		} else if (object@rescue.ext <= 0L){
			return("extension length must be a positive integer")	
		}
	}
		
	if (length(object@dedup)!=1L || !is.logical(object@dedup)) { 
		return("duplicate removal specification must be a logical scalar")
	}
	if (length(object@forward)!=1L || !is.logical(object@forward)) { 
		return("forward strand specification must be a logical scalar")
	}
	if (length(object@minq)!=1L || !is.numeric(object@minq)) { 
		return("minimum mapping quality must be a numeric scalar")
	}
	return(TRUE)
})

setMethod("initialize", signature("readParam"), function(.Object, ...) {
	value <- callNextMethod()
	validObject(value)
	value
})

setMethod("$", signature("readParam"), function(x, name) { 
	slot(x, name)
})

setMethod("show", signature("readParam"), function(object) {
	cat("    ", switch(object@pe,
 	   none="Extracting reads in single-end mode",
	   both="Extracting reads in paired-end mode",
	   first="Extracting the first read of each pair",
	   second="Extracting the second read of each pair"), "\n", sep="")

	if (object@pe=="both") { 
		cat("        Maximum allowed distance between paired reads is", object@max.frag, "bp\n")
		cat("        Rescuing of improperly paired reads is", 
			ifelse(object@rescue.pairs, "enabled", "disabled"), "\n")
		if (object@rescue.pairs) {
			cat("            Extension length for rescued reads is", object@rescue.ext, "bp\n")
		}
	}

	cat("    Duplicate removal is turned", ifelse(object@dedup, "on", "off"), "\n")
	if (is.na(object@minq)) { 
		cat("    No minimum threshold is set on the mapping score\n")
	} else {
		cat("    Minimum allowed mapping score is", object@minq, "\n")
	}

	if (is.na(object@forward)) { 
		cat("    Strand specificity is turned off\n")
	} else {
		cat("    Reads are extracted from the", ifelse(object@forward, "forward", "reverse"), "strand only\n")
	}

	rl <- length(object@restrict)
	if (rl) { 
		cat("    Read extraction is limited to", rl, ifelse(rl==1L, "sequence\n", "sequences\n"))
	} else {
		cat("    No restrictions are placed on read extraction\n")
	}

	dl <- length(object@discard)
	if (dl) { 
		cat("    Reads in", dl, ifelse(dl==1L, "region", "regions"), "will be discarded\n")
	} else {
		cat("    No regions are specified to discard reads\n")
	}
})

readParam <- function(pe="none", max.frag=500, rescue.pairs=FALSE,
	rescue.ext=NA, dedup=FALSE, minq=NA, forward=NA, restrict=NULL, discard=GRanges())
# This creates a SimpleList of parameter objects, specifying
# how reads should be extracted from the BAM files. The aim is
# to synchronize read loading throughout the package, such that
# you don't have to manually respecify them in each function.
#
# written by Aaron Lun
# created 1 September 2014
{
	max.frag <- as.integer(max.frag)
	rescue.pairs <- as.logical(rescue.pairs)
	if (rescue.pairs) {
		if (is.na(rescue.ext)) { stop("need to specify extension length for rescued reads") }
 	}	
	rescue.ext <- as.integer(rescue.ext)

	dedup <- as.logical(dedup)
	forward <- as.logical(forward)
	minq <- as.integer(minq)
	restrict <- as.character(restrict) 
	new("readParam", pe=pe, max.frag=max.frag, rescue.pairs=rescue.pairs,
		rescue.ext=rescue.ext, dedup=dedup, forward=forward, minq=minq, 
		restrict=restrict, discard=discard)
}

setGeneric("reform", function(x, ...) { standardGeneric("reform") })
setMethod("reform", signature("readParam"), function(x, ...) {
	incoming <- list(...)
	sn <- slotNames(x)
	for (sx in names(incoming)) {
		val <- incoming[[sx]]
		sx <- match.arg(sx, sn)
		incoming[[sx]] <- switch(sx, 
			max.frag=as.integer(val),
			rescue.pairs=as.logical(val),
			rescue.ext=as.integer(val),
			dedup=as.logical(val),
			forward=as.logical(val),
			minq=as.integer(val),
			restrict=as.character(val),
			val)
	}
	do.call(initialize, c(x, incoming))
})

