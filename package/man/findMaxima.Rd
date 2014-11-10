\name{findMaxima}
\alias{findMaxima}

\title{Find local maxima}
\description{Find the local maxima for a given set of genomic regions.}

\usage{
findMaxima(regions, range, metric)
}

\arguments{
\item{regions}{a \code{GRanges} object}
\item{range}{an integer scalar specifying the range of surrounding regions to consider as local}
\item{metric}{a numeric vector of values for which the local maxima is found, defaults to the average count}
}

\details{
For each region in \code{regions}, this function will examine all regions within \code{range} on either side.
It will then determine if the current region has the maximum value of \code{metric} across this range.
A typical metric to maximize might be the sum of counts or the average abundance across all libraries.

Preferably, \code{regions} should contain regularly sized and spaced windows or bins, e.g., from \code{\link{windowCounts}}.
The sensibility of using this function for arbitrary regions is left to the user.
In particular, the algorithm will not support nested regions and will fail correspondingly if any are detected.
}

\value{
A logical vector indicating whether each region in \code{regions} is a local maxima.
}

\seealso{
\code{\link{windowCounts}},
\code{\link{aveLogCPM}}
}

\author{
Aaron Lun
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, filter=1)
regions <- rowData(data)
metric <- edgeR::aveLogCPM(asDGEList(data))
findMaxima(regions, range=10, metric=metric)
findMaxima(regions, range=50, metric=metric)
findMaxima(regions, range=100, metric=metric)

findMaxima(regions, range=10, metric=runif(length(regions)))
findMaxima(regions, range=50, metric=runif(length(regions)))
findMaxima(regions, range=100, metric=runif(length(regions)))
}

\keyword{filtering}
