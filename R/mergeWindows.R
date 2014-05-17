mergeWindows <- function(regions, tol, sign=NULL, max.width=NULL) 
# This function merges the adjacent windows if they lie within 'tol' of each other,
# Any abundance filtering should be done beforehand. Negative values of tol refer
# to a minimum overlap. A value of zero means that the windows must be adjacent
# (i.e. non-overlapping and contiguous).
# 
# written by Aaron Lun
# 30 July 2013
{
	tol<-as.integer(tol+0.5)
	max.width<-as.integer(max.width+0.5)
	o<-GenomicRanges::order(regions)
	regions<-regions[o]
	if (is.null(sign)) { 
		sign<-logical(length(regions)) 
	} else {
		sign<-sign[o]
	}

	# Running the merge.
	out<-.Call("R_merge", as.integer(seqnames(regions)), start(regions), end(regions), sign, 
			tol, max.width, PACKAGE="csaw")
	if (is.character(out)) { stop(out) }
	
	# Reporting. Indices correspond with positions in 'clustered'.
	out[[1]][o]<-out[[1]]
	clustered<-GRanges(levels(seqnames(regions))[out[[2]]], IRanges(out[[3]], out[[4]]),
			seqinfo=Seqinfo(seqlevels(regions), seqlengths(regions)))
	return(list(id=out[[1]], region=clustered))
}
