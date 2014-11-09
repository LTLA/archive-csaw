findMaxima <- function(regions, range, metric)
# This function finds the maximum window in 'data', given a range
# around which the maxima is to be considered. The 'metric' is,
# by default, the average count, but others can be used if supplied.
#
# written by Aaron Lun
# Created 9 November 2014.
{
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
