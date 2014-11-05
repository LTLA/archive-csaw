\name{scaledAverage}
\alias{scaledAverage}

\title{Scaled average abundance}
\description{Compute the scaled average abundance for each feature.}

\usage{
scaledAverage(y, scale=1, prior.count=NULL, ...)
}

\arguments{
\item{y}{a \code{DGEList} object}
\item{scale}{a numeric vector indicating the magnitude with which each abundance is to be downscaled}
\item{prior.count}{a numeric scalar specifying the prior count to add}
\item{...}{other arguments, to be passed to \code{\link{aveLogCPM}}}
}

\details{
This function computes the average abundance of each feature in \code{y}, and downscales it according to \code{scale}.
For example, if \code{scale=2}, the average count is halved, i.e., the returned abundances are decreased by 1 (as they are log2-transformed value).
The aim is to set \code{scale} based on the relative width of regions, to allow abundances to be compared between regions of differing widths.
Widths can be obtained using \code{\link{getWidths}}.

Some subtlety is necessary regarding the treatment of the \code{prior.count}.
Specifically, the prior used in \code{\link{aveLogCPM}} is automatically increased when \code{scale} is larger.
This ensures that the effective prior is the same after the abundance is scaled down.
Otherwise, the use of the same prior would incorrectly result in a smaller abundance for larger regions.

Note that the adjustment for width assumes that reads are uniformly distributed throughout each region.
This is reasonable for most background regions, but may not be for enriched regions.
When the distribution is highly heterogeneous, the downscaled abundance of a large region will not be an accurate representation of the abundance of the smaller regions nested within.

For consistency, the \code{prior.count} is set to the default value of \code{\link{aveLogCPM.DGEList}}, if it is not otherwise specified.
If a non-default value is used, make sure that it is the same for all calls to \code{scaledAverage}.
This ensures that comparisons between the returned values are valid.
}

\value{
A numeric vector of scaled abundances, with one entry for each row of \code{y}.
}

\seealso{
\code{\link{getWidths}},
\code{\link{aveLogCPM}}
}

\author{
Aaron Lun
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
size1 <- 50
data1 <- windowCounts(bamFiles, width=size1, filter=1)
size2 <- 500
data2 <- windowCounts(bamFiles, width=size2, filter=1)

# Adjusting by `scale`, based on expected sizes.
scaledAverage(asDGEList(data1), scale=1)
scaledAverage(asDGEList(data2), scale=size2/size1)

# Alternatively, can use observed sizes.
relative <- getWidths(data2)/getWidths(data1)
scaledAverage(asDGEList(data1))
scaledAverage(asDGEList(data2), scale=relative)

# Need to make sure the same prior is used, if non-default.
pc <- 5
scaledAverage(asDGEList(data1), prior.count=pc)
scaledAverage(asDGEList(data2), scale=relative, prior.count=pc)
}
