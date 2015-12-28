csawUsersGuide <- function(view=TRUE)
# Find, and optionally show, the csaw user's guide 
# 
# written by Aaron Lun, based on equivalent from edgeR by Gordon Smyth
# 21 June 2014.
{
	ugloc <- system.file("doc", "csawUserGuide.pdf", package="csaw")
	if (view) {
		if(.Platform$OS.type == "windows") {
			shell.exec(ugloc)
		} else {
			system(paste(Sys.getenv("R_PDFVIEWER"),ugloc,"&"))
		}
	}
	return(ugloc)
}
