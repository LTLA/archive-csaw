findMaxima <- function(regions, range, metric, ignore.strand=TRUE)
# This function finds the maximum window in 'data', given a range
# around which the maxima is to be considered. The 'metric' is,
# by default, the average count, but others can be used if supplied.
#
# written by Aaron Lun
# created 9 November 2014.
# last modified 10 February 2015.
{
	strs <- strand(regions)
	if (!ignore.strand && length(runValue(strs))!=1) {
		# Strand-specific maxima identification.
		forward <- as.logical(strs=="+")
		reverse <- as.logical(strs=="-")
		neither <- as.logical(strs=="*")
		out <- logical(length(regions))
		if (any(forward)) { out[forward] <- Recall(regions=regions[forward], range=range, metric=metric[forward], ignore.strand=TRUE) }
		if (any(reverse)) { out[reverse] <- Recall(regions=regions[reverse], range=range, metric=metric[reverse], ignore.strand=TRUE) }
		if (any(neither)) { out[neither] <- Recall(regions=regions[neither], range=range, metric=metric[neither], ignore.strand=TRUE) }
		return(out)
	}

	chrs <- as.integer(seqnames(regions))
	starts <- start(regions)
	ends <- end(regions)
	o <- order(chrs, starts, ends)

	if (length(metric)!=length(regions)) { stop("one metric must be supplied per region") }
	if (!is.double(metric)) { metric <- as.double(metric) }
	stopifnot(!any(is.na(metric)))
	if (!is.integer(range)) { range <- as.integer(range) }

	out <- .Call(cxx_find_maxima, chrs[o], starts[o], ends[o], metric[o], range)
	if (is.character(out)) { stop(out) }
	out[o] <- out
	return(out)
}
