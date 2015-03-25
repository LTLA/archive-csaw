\name{overlapStats}
\alias{combineOverlaps}
\alias{getBestOverlaps}
\alias{summitOverlaps}

\title{Compute overlap statistics}
\description{Compute assorted statistics for overlaps between windows and regions in a Hits object.}

\usage{
combineOverlaps(olap, tab, relation.weights=NULL, test.weights=NULL, ...)
getBestOverlaps(olap, tab, relation.weights=NULL, test.weights=NULL, ...)
summitOverlaps(olap, region.best, test.summits)
}

\arguments{
\item{olap}{a \code{Hits} object produced by \code{\link{findOverlaps}}, containing overlaps between regions (query) and windows (subject)}
\item{tab}{a dataframe of DE results for each window}
\item{relation.weights}{a numeric vector specifying weights for each overlapped window}
\item{test.weights}{a numeric vector specifying weights for each individual window}
\item{...}{other arguments to be passed to the wrapped functions}
\item{region.best}{an integer vector specifying the window index that is the summit for each region}
\item{test.summits}{a logical vector specifying whether an individual window is a summit, or an integer vector indicating the windows that correspond to summits}
}

\details{
These functions provide convenient wrappers around \code{\link{combineTests}}, \code{\link{getBestTest}} and \code{\link{upweightSummit}}.
They accept \code{Hits} objects produced by running \code{\link{findOverlaps}} between windows and some pre-specified regions.
Each set of windows overlapping a region is defined as a cluster to compute the various statistics.

A wrapper is necessary as a window may overlap multiple regions.
If so, the multiple instances of that window are defined as distinct ``overlapped'' windows, where each overlapped window is assigned to a different region.
In contrast, the ``individual'' window just refers to the window itself, regardless of what it overlaps.

The \code{relation.weights} argument refers to the weights for each region-window relationship, i.e., each row of \code{olap}.
This allows for different weights to be assigned to the same window in different regions.
The \code{test.weights} argument is the weight of the window itself, and is the same regardless of the region.
If both are specified, \code{relation.weights} takes precedence.

For \code{summitOverlaps}, the \code{region.best} argument is designed to accept the \code{best} field in the output of \code{getBestOverlaps} (run with \code{by.pval=FALSE}).
This indicates whether an overlapped window is a summit within each region.
In contrast, the \code{test.weights} argument indicates whether an individual window is a summit, e.g., from \code{\link{findMaxima}}.
}

\value{
For \code{combineOverlaps} and \code{getBestOverlaps}, a dataframe is returned from their respective wrapped functions.
Each row of the dataframe corresponds to a region, where regions without overlapped windows are assigned \code{NA} values.

For \code{summitOverlaps}, a numeric vector of weights is produced.
This can be used as \code{relation.weights} in the other two functions.
}

\seealso{
\code{\link{combineTests}},
\code{\link{getBestTest}},
\code{\link{upweightSummit}}
}

\author{
Aaron Lun
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=1, filter=1)
of.interest <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(200, 1000, 700, 1500)))

# Making some mock results.
N <- nrow(data)
mock <- data.frame(logFC=rnorm(N), PValue=runif(N), logCPM=rnorm(N))

olap <- findOverlaps(of.interest, rowRanges(data))
combineOverlaps(olap, mock)
getBestOverlaps(olap, mock)

# See what happens when you don't get many overlaps.
getBestOverlaps(olap[1,], mock)
combineOverlaps(olap[2,], mock)

# Weighting example, with window-specific weights.
window.weights <- runif(N) 
comb <- combineOverlaps(olap, mock, test.weights=window.weights)

# Weighting example, with relation-specific weights.
best.by.ave <- getBestOverlaps(olap, mock, by.pval=FALSE)
w <- summitOverlaps(olap, region.best=best.by.ave$best)
head(w)
stopifnot(length(w)==length(olap))
combineOverlaps(olap, mock, relation.weights=w)

# Running summitOverlaps for window-specific summits
# (output is still relation-specific weights, though).
is.summit <- findMaxima(rowRanges(data), range=100, metric=mock$logCPM)
w <- summitOverlaps(olap, test.summits=is.summit)
head(w)
}

