maximizeCcf <- function(profile, ignore=100) 
# A quick and dirty function to maximize the CCF to get the 
# average fragment length, while ignoring the phantom peak.
#
# written by Aaron Lun
# created 15 August 2015
{
    ignore <- as.integer(ignore)
    if (ignore >= length(profile)) { stop("cannot ignore the entire profile") }
    if (ignore < 0L) { stop("ignore length must be positive") }
    which.max(profile[(ignore+1):length(profile)]) + ignore - 1L
}
