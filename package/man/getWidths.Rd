\name{getWidths}
\alias{getWidths}

\title{Get region widths}
\description{Get the widths of the read counting interval for each region.}

\usage{
getWidths(data)
}

\arguments{
\item{data}{a \code{SummarizedExperiment} object, produced by \code{\link{windowCounts}} or \code{\link{regionCounts}}}
}

\details{
Widths of all regions are increased by the average fragment length during the calculations.
This is because each count represents the number of (imputed) fragments overlapping each region.
Thus, a 1 bp window has an effective width that includes the average length of each fragment. 
The fragment length is taken from \code{data$final.ext} for single-end data.

% It's a direct addition of width and fragment length, minus 1. If you have a 100 bp interval,
% and you slide a 50 bp fragment along the genome until it hits the interval, there
% are 149 positions where the fragment overlaps the interval. So, the effective width
% is actually 149 bp. This logic holds regardless of which strand the read was from, so
% there's no need to adjust things if data$param$forward is not NA.

For paired-end data, the average fragment length should be the median of the values obtained with \code{\link{getPESizes}}.
This can be specified in the \code{rescue.ext} slot of each \code{data$param}.
Alternatively, if rescuing is disabled, the length will be extracted from the value of \code{data$ext}.
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

bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, ext=c(50, 100), filter=1)
getWidths(data)

bamFile <- system.file("exdata", "pet.bam", package="csaw")
data <- windowCounts(bamFile, param=readParam(pe="both"), filter=1)
getWidths(data)

# Avoid warning by defining 'rescue.ext'.
exptData(data)$param <- reformList(exptData(data)$param, rescue.ext=200)
getWidths(data)
}

\keyword{filtering}
