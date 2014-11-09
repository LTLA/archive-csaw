\name{findMaxima}
\alias{findMaxima}

\title{Find local maxima}
\description{Find the local maxima for a given set of genomic regions.}

\usage{
findMaxima(data, range, metric=NULL)
}

\arguments{
\item{data}{a \code{SummarizedExperiment} object, produced by \code{\link{windowCounts}}}
\item{range}{an integer scalar specifying the range of surrounding regions to consider as local}
\item{metric}{a numeric vector of values for which the local maxima is found, defaults to the average count}
}

\details{
For each region in \code{data}, this function will examine all regions within \code{range} on either side.
It will then determine if the current region has the maximum value of \code{metric} across this range.
By default, \code{metric} is set to the average abundance (as computed by \code{\link{aveLogCPM}}) if not otherwise specified.

Preferably, \code{data} should be constructed using regularly sized and spaced windows or bins, e.g., from \code{\link{windowCounts}}.
The sensibility of using this function for arbitrary regions is left to the user.
The algorithm will not support nested regions and will fail correspondingly if any are detected.
}

\value{
A logical vector indicating whether each region in \code{data} is a local maxima.
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
findMaxima(data, range=10)
findMaxima(data, range=50)
findMaxima(data, range=100)

findMaxima(data, range=10, metric=runif(nrow(data)))
findMaxima(data, range=50, metric=runif(nrow(data)))
findMaxima(data, range=100, metric=runif(nrow(data)))
}
