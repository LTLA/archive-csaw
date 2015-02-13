makeExtVector <- function(ext, final.ext=NULL) 
# This is a convenience function that manufactures the extension vector, with
# the final extension length stuffed into the attributes as 'final.ext'.
#
# written by Aaron Lun
# created 13 February 2015
{
	if (is.null(final.ext)) { final.ext <- mean(ext) }
	final.ext <- as.integer(final.ext)
	ext <- as.integer(ext)
	attributes(ext)$final.ext <- final.ext
	class(ext) <- "extVector"
	return(ext)
}

`[.extVector` <- function(x, i, ...) {
	attrs <- attributes(x)
	out <- unclass(x)
	out <- out[i]
	attributes(out) <- attrs
	out
}
