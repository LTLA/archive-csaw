\name{mergeWindows}
\alias{mergeWindows}

\title{Merge windows into clusters}

\description{Uses a simple single-linkage approach to merge adjacent or overlapping windows into clusters.}

\usage{
mergeWindows(regions, tol, sign=NULL, max.width=NULL) 
}

\arguments{
	\item{regions}{a \code{GRanges} object}
	\item{tol}{a numeric scalar specifying the maximum distance between adjacent windows}
	\item{sign}{a logical vector specifying whether each window has a positive log-FC}
	\item{max.width}{a numeric scalar specifying the maximum size of merged intervals}
}

\value{
	A list containing \code{id}, an integer vector containing the cluster ID for each window; and \code{region},
	  a \code{GRanges} object containing the start/stop coordinates for each cluster of windows.
}

\details{
Windows are merged if the gap between the end of one window and the start of
the next is no greater than \code{tol}. Adjacent windows can then be chained
together to build a cluster of windows across the linear genome. A value of
zero means that the windows must be contiguous whereas negative values specify
minimum overlaps.

If \code{sign!=NULL}, windows are only merged if they have the same sign of the
log-FC and are not separated by windows with opposite log-FC values. This can
be useful when summarizing adjacent regions. However, it is not recommended for
routine clustering in differential analyses as the resulting clusters will not
be independent of the p-value.

The \code{max.width} parameter can be specified to avoid the formation of
excessively large clusters when many adjacent regions are present. For typical
ChIP-seq datasets, suggested values range from 2000 to 10000 bp. Clusters are
split if they exceed this size, albeit in a rather \emph{ad hoc} manner. If
\code{NULL}, no limits are placed on the maximum size.

The tolerance should reflect the minimum distance at which two regions of
enrichment are considered separate. If two windows are more than \code{tol}
apart, they \emph{will} be placed into separate clusters. In contrast, the
\code{max.width} value reflects the maximum distance at which two windows can be
considered part of the same region. 

Note that in the output, the cluster ID reported in \code{id} corresponds to
the index of the cluster coordinates in the input \code{region}.
}

\examples{
x<-round(runif(10, 100, 1000))
gr<-GRanges(rep("chrA", 10), IRanges(x, x+40))
mergeWindows(gr, 1)
mergeWindows(gr, 10)
mergeWindows(gr, 100)
mergeWindows(gr, 100, sign=rep(c(TRUE, FALSE), 5))
}

\seealso{
\code{\link{combineTests}}, \code{\link{windowCounts}}
}

\author{Aaron Lun}

\keyword{clustering}
