# Assorted convenience functions to handle and manipulate lists of readParam
# objects.
#
# written by Aaron Lun
# created 12 December 2014
# last modified 7 February 2015

reformList <- function(paramlist, ...) {
	lapply(paramlist, FUN=reform, ...)
}

checkList <- function(paramlist) {
	if (length(paramlist)>1L) { 
		different <- list(character(0))
		sn <- slotNames(paramlist[[1]])
		for (i in 2:length(paramlist)) { 
			for (x in sn) {
				if (!identical(	slot(paramlist[[i]], x), slot(paramlist[[1]], x))) {
					different[[length(different)+1L]] <- x
				}
			}
		}
		return(unlist(different))
	}
	return(character(0))
}
