findMaxima <- function(data, range, metric=NULL)
# This function finds the maximum window in 'data', given a range
# around which the maxima is to be considered. The 'metric' is,
# by default, the average count, but others can be used if supplied.
#
# written by Aaron Lun
# Created 9 November 2014.
{
	regions <- rowData(data)
	chrs <- as.integer(seqnames(regions))
	starts <- start(regions)
	ends <- end(regions)
	o <- order(chrs, starts, ends)

	if (is.null(metric)) { metric <- aveLogCPM(asDGEList(data)) }
	if (!is.double(metric)) { metric <- as.double(metric) }
	stopifnot(!any(is.na(metric)))
	if (!is.integer(range)) { range <- as.integer(range) }

	out <- .Call(cxx_find_maxima, chrs[o], starts[o], ends[o], metric[o], range)
	if (is.character(out)) { stop(out) }
	out[o] <- out
	return(out)
}
