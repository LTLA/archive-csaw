\name{getWidths}
\alias{getWidths}

\title{Get region widths}
\description{Get the widths of the read counting interval for each region.}

\usage{
getWidths(data)
}

\arguments{
\item{data}{a \code{RangedSummarizedExperiment} object, produced by \code{\link{windowCounts}} or \code{\link{regionCounts}}}
}

\details{
Widths of all regions are increased by the average fragment length during the calculations.
This is because each count represents the number of (imputed) fragments overlapping each region.
Thus, a 1 bp window has an effective width that includes the average length of each fragment. 

% It's a direct addition of width and fragment length, minus 1. If you have a 100 bp interval,
% and you slide a 50 bp fragment along the genome until it hits the interval, there
% are 149 positions where the fragment overlaps the interval. So, the effective width
% is actually 149 bp (compared to the case where ext=1; imagine the interval spanned by the
% midpoints of the fragments, or by the start or end of the fragment for forward or reverse
% strand reads, respectively). This logic holds regardless of which strand the read was from, 
% so there's no need to adjust things if data$param$forward is not NA.

The fragment length is taken from \code{metadata(data)$final.ext}, if it is not \code{NA}.
Otherwise, it is taken from \code{data$ext}.
If the fragment lengths are different between libraries, the average is used to computed the effective width of the window.
For paired-end data, \code{data$ext} should be an average of the inferred fragment sizes, e.g., obtained with \code{\link{getPESizes}}.

If \code{final.ext} is \code{NA} and any of \code{ext} are \code{NA}, the function will try to extract the read lengths in \code{data$rlen}.
This is because \code{NA} values of \code{ext} correspond to the use of unextended reads in \code{\link{windowCounts}} and \code{\link{regionCounts}}.
The likely read lengths are automatically computed in each function but can also be set manually.
}

\value{
An integer vector containing the effective width, in base pairs, of each region.
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}}
}

\author{
Aaron Lun
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, filter=1)
getWidths(data)

data <- windowCounts(bamFiles, ext=c(50, 100), filter=1)
getWidths(data)

data <- windowCounts(bamFiles, ext=DataFrame(c(50, 100), NA), filter=1)
getWidths(data)

# Using the automatically computed 'rlen'.
data <- windowCounts(bamFiles, ext=c(NA, 100), filter=1)
getWidths(data)
data$rlen <- 200 # Manually defining it, if required.
getWidths(data)

# Paired-end data also takes the fragment length from 'ext'.
bamFile <- system.file("exdata", "pet.bam", package="csaw")
data <- windowCounts(bamFile, param=readParam(pe="both"), filter=1)
getWidths(data)
data$ext <- 200 # Again, manual definition is accepted.
getWidths(data)
}

\keyword{filtering}
