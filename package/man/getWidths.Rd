\name{getWidths}
\alias{getWidths}

\title{Get region widths}
\description{Get the widths of the read counting interval for each region.}

\usage{
getWidths(data, len=NULL)
}

\arguments{
\item{data}{a \code{SummarizedExperiment} object, produced by \code{\link{windowCounts}} or \code{\link{regionCounts}}}
\item{len}{an integer scalar of vector of fragment lengths for each library}
}

\details{
Widths of all regions are increased by the average fragment length during the calculations.
This is because each count represents the number of (imputed) fragments overlapping each region.
Thus, a 1 bp window has an effective width that includes the average length of each fragment. 

% It's a direct addition of width and fragment length, minus 1. If you have a 100 bp interval,
% and you slide a 50 bp fragment along the genome until it hits the interval, there
% are 149 positions where the fragment overlaps the interval. So, the effective width
% is actually 149 bp. This logic holds regardless of which strand the read was from, so
% there's no need to adjust things if data$param$forward is not NA.

The fragment length is taken from \code{exptData(data)$final.ext}, if it is not \code{NA}.
Otherwise, for single-end data, it is taken from \code{data$ext}.
For paired-end data, the average fragment length should be the median of the values obtained with \code{\link{getPESizes}}.
This is taken from the \code{rescue.ext} slot of each \code{data$param}.

If the fragment lengths are different between libraries, the average is used to computed the effective width of the window.
If \code{final.ext} is \code{NA} and any of \code{ext} or \code{rescue.ext} are \code{NA} for the relevant libraries, the function will throw an error.
Users should supply valid fragment lengths in \code{len} to avoid this.
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

# Avoid error by defining 'len', or by rescuing with 'rescue.ext'.
bamFile <- system.file("exdata", "pet.bam", package="csaw")
data <- windowCounts(bamFile, param=readParam(pe="both"), filter=1)
try(getWidths(data))

getWidths(data, len=200)
data <- windowCounts(bamFile, param=readParam(pe="both", rescue.ext=200), filter=1)
getWidths(data)
}

\keyword{filtering}
