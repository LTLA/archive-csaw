# This defines the readParam class, which specifies the universal
# parameters for read loading in csaw. A class is used here so that
# there are continuous validity checks on the list values.

setClass("readParam", representation(pet="character", 
	max.frag="integer", rescue.pairs="logical", rescue.ext="integer",
	dedup="logical", minq="integer", 
	restrict="character", discard="GRanges"))

setValidity("readParam", function(object) {
    if (length(object@pet)!=1L || ! object@pet %in%	c("none", "both", "first", "second")) { 
		return("PET specification must be a character scalar of 'none', 'both', 'first' or 'second'") 
	}
   	if (length(object@max.frag)!=1L || object@max.frag <= 0L) {
		return("maximum fragment specifier must be a positive integer")
	} 
	if (length(object@rescue.pairs)!=1L) {
		return("pair rescue specification must be a logical scalar")
	}
	if (length(object@rescue.ext)!=1L || object@rescue.ext <= 0L){
		return("rescue extension length must be a positive integer")	
	}
	if (length(object@dedup)!=1L || !is.logical(object@dedup)) { 
		return("duplicate removal specification must be a logical scalar")
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

setMethod("$<-", signature("readParam"), function(x, name, value) {
	slot(x, name) <- value
	validObject(x)
	x	
})

setMethod("show", signature("readParam"), function(object) {
	cat(switch(object@pet,
 	   none="Extracting reads in single-end mode",
	   both="Extracting reads in paired-end mode",
	   first="Extracting the first read of each pair",
	   second="Extracting the second read of each pair"), "\n")
	if (object@pet=="both") { 
		cat("  Maximum allowed distance between paired reads is", object@max.frag, "\n")
		cat("  Rescuing of improperly paired reads is", 
			ifelse(object@rescue.pairs, "enabled", "disabled"), "\n")
		if (object@rescue.pairs) {
			cat("    Extension length for rescued reads is", object@rescue.ext, "\n")
		}
	}
	cat("Duplicate removal is turned", ifelse(object@dedup, "on", "off"), "\n")
	if (is.na(object@minq)) { 
		cat("No minimum threshold is set on the mapping score\n")
	} else {
		cat("Minimum allowed mapping score is", object@minq, "\n")
	}
	if (length(object@restrict)) { 
		cat("Read extraction is limited to", length(object@restrict), "sequences\n")
	} else {
		cat("No restrictions are placed on read extraction\n")
	}
	if (length(object@discard)) { 
		cat("Reads in", length(object@discard), "regions will be discarded\n")
	} else {
		cat("No regions are specified to discard reads\n")
	}
})

readParam <- function(pet=c("none", "both", "first", "second"), max.frag=500, rescue.pairs=FALSE,
	rescue.ext=200, dedup=FALSE, minq=NA, restrict=NULL, discard=GRanges())
# This creates a SimpleList of parameter objects, specifying
# how reads should be extracted from the BAM files. The aim is
# to synchronize read loading throughout the package, such that
# you don't have to manually respecify them in each function.
#
# written by Aaron Lun
# 1 September 2014
{
	pet <- match.arg(pet)
	max.frag <- as.integer(max.frag)
	rescue.pairs <- as.logical(rescue.pairs)
	rescue.ext <- as.integer(rescue.ext)
	dedup <- as.logical(dedup)
	minq <- as.integer(minq)
	restrict <- as.character(restrict) 
	new("readParam", pet=pet, max.frag=max.frag, rescue.pairs=rescue.pairs,
		rescue.ext=rescue.ext, dedup=dedup, minq=minq, restrict=restrict, 
		discard=discard)
}

