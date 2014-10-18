apkb <- function(data, prior.count=10, pet.len=NULL, widths=NULL, ...) 
# This calculates the average abundance-per-kilobase. It basically represents
# the average log-CPM, adjusted for the window sizes of every region involved.
# This requires a careful touch because of the treatment of the prior count;
# you need to supply a large prior count if you're downscaling later.
# 
# written by Aaron Lun 
# 18 October 2014
{
	if (is.null(widths)) { 
		is.pet <- exptData(data)$param$pet=="both"
		if (is.pet) {
			if (is.null(pet.len)) { stop("pet.len must be specified for paired-end data") }
			pet.len <- as.integer(pet.len)
			if (pet.len <= 0L) { stop("pet.len must be a positive integer") }
			frag.len <- pet.len
		} else {
			frag.len <- exptData(data)$ext
		}
		widths <- width(rowData(data)) + frag.len - 1L
	}

	# Checking the supplied values.
   	if (length(widths)!=nrow(data)) { 
		stop("length of widths should be equal to elements in data")
	}	
	if (any(widths <= 0L)) { 
		warning("computed widths should be positive integers")
		widths <- pmax(0L, widths)
	}

	# Computing abundances.
	adjustment <- widths/1000
	new.prior <- adjustment * prior.count
	ab <- aveLogCPM(assay(data), lib.size=data$totals, prior.count=new.prior, ...)	
	ab - log2(adjustment)
}
