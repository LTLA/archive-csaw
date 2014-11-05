scaledAverage <- function(y, scale=1, prior.count=NULL, ...)
# This computes the scaled average abundance, with some finesse to deal with
# the interaction between scaling and the prior count. The `scale` factor
# represents the downscaling factor for the abundances, so the prior count
# has to be upscaled during the calculation itself.
#
# written by Aaron Lun
# 5 November 2014
{
	if (is.null(prior.count)) { prior.count <- formals(aveLogCPM.DGEList)$prior.count }
	aveLogCPM(y, prior.count=scale*prior.count, ...) - log2(scale)
}

getWidths <- function(data, pet.len=NULL) 
# This computes the effective width of the data in the SummarizedExperiment
# object. This is done by accounting for the effect of read extension; or,
# for paired end data, the median fragment length.
#
# written by Aaron Lun
# 5 November 2014
{
	is.pet <- exptData(data)$param$pet=="both"
	if (is.pet) {
		if (is.null(pet.len)) { stop("pet.len must be specified for paired-end data") }
		pet.len <- as.integer(pet.len)
		if (pet.len <= 0L) { stop("pet.len must be a positive integer") }
		frag.len <- pet.len
	} else {
		frag.len <- exptData(data)$ext
	}
	width(rowData(data)) + frag.len - 1L
}

