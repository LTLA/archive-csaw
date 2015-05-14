# This defines the readParam class, which specifies the universal
# parameters for read loading in csaw. A class is used here so that
# there are continuous validity checks on the list values.

setClass("readParam", representation(pe="character", 
	max.frag="integer", rescue.ext="integer", fast.pe="logical",
	dedup="logical", minq="integer", forward="logical", 
	restrict="character", discard="GRanges"))

setValidity("readParam", function(object) {
    if (length(object@pe)!=1L || ! object@pe %in%	c("none", "both", "first", "second")) { 
		return("PE specification must be a character scalar of 'none', 'both', 'first' or 'second'") 
	}
   	if (length(object@max.frag)!=1L || object@max.frag <= 0L) {
		return("maximum fragment specifier must be a positive integer")
	} 

	if (length(object@rescue.ext)!=1L) {
		return("rescue specification must be a integer scalar")
	} else if (!is.na(object@rescue.ext) && object@rescue.ext <= 0L){
		return("extension length must be a positive integer")	
	}
	if (length(object@fast.pe)!=1L || !is.logical(object@fast.pe)) { 
		return("fast PE extraction flag must be a logical scalar")
	}
		
	if (length(object@dedup)!=1L || !is.logical(object@dedup)) { 
		return("duplicate removal specification must be a logical scalar")
	}
	if (length(object@minq)!=1L || !is.numeric(object@minq)) { 
		return("minimum mapping quality must be a numeric scalar")
	}

	if (length(object@forward)>1L || !is.logical(object@forward)) { 
		return("forward strand specification must be logical")
	} else if ((length(object@forward)==0L || !is.na(object@forward)) && object@pe == "both") {
		stop("strand-specific extraction makes no sense for paired-end data")
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
   		if (object@fast.pe) {
			cat("        Fast PE extraction enabled, ignoring other settings\n")
			return(invisible(NULL))
		}	
		if (!is.na(object@rescue.ext)) {
			cat("        Rescuing of improperly paired reads is enabled\n")
			cat("            Extension length for rescued reads is", object@rescue.ext, "bp\n")
		} else {
			cat("        Rescuing of improperly paired reads is disabled\n")
		}
	}

	cat("    Duplicate removal is turned", ifelse(object@dedup, "on", "off"), "\n")
	if (is.na(object@minq)) { 
		cat("    No minimum threshold is set on the mapping score\n")
	} else {
		cat("    Minimum allowed mapping score is", object@minq, "\n")
	}

	if (length(object@forward)==0L) { 
		cat("    Reads are extracted from either strand separately\n")
	} else if (is.na(object@forward)) { 
		cat("    Reads are extracted from both strands\n")
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

readParam <- function(pe="none", max.frag=500, rescue.ext=NA, fast.pe=FALSE, 
	dedup=FALSE, minq=NA, forward=NA, restrict=NULL, discard=GRanges())
# This creates a list of parameters, formally represented as a readParam
# object, specifying how reads should be extracted from the BAM files. The
# aim is to synchronize read loading throughout the package, such that
# you don't have to manually respecify them in each function.
#
# written by Aaron Lun
# created 1 September 2014
{
	max.frag <- as.integer(max.frag)
	rescue.ext <- as.integer(rescue.ext)
	fast.pe <- as.logical(fast.pe)
	dedup <- as.logical(dedup)
	forward <- as.logical(forward)
	minq <- as.integer(minq)
	restrict <- as.character(restrict) 
	new("readParam", pe=pe, max.frag=max.frag, fast.pe=fast.pe,
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
			rescue.ext=as.integer(val),
			fast.pe=as.logical(val),
			dedup=as.logical(val),
			forward=as.logical(val),
			minq=as.integer(val),
			restrict=as.character(val),
			val)
	}
	do.call(initialize, c(x, incoming))
})

