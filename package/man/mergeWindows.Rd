\name{mergeWindows}
\alias{mergeWindows}

\title{Merge windows into clusters}

\description{Uses a simple single-linkage approach to merge adjacent or overlapping windows into clusters.}

\usage{
mergeWindows(regions, tol, sign=NULL, max.width=NULL, ignore.strand=TRUE) 
}

\arguments{
	\item{regions}{a \code{GRanges} object}
	\item{tol}{a numeric scalar specifying the maximum distance between adjacent windows}
	\item{sign}{a logical vector specifying whether each window has a positive log-FC}
	\item{max.width}{a numeric scalar specifying the maximum size of merged intervals}
	\item{ignore.strand}{a logical scalar indicating whether to consider the strandedness of \code{regions}}
}

\value{
	A list containing \code{id}, an integer vector containing the cluster ID for each window; and \code{region},
	  a \code{GRanges} object containing the start/stop coordinates for each cluster of windows.
}

\details{
Windows are merged if the gap between the end of one window and the start of the next is no greater than \code{tol}. 
Adjacent windows can then be chained together to build a cluster of windows across the linear genome. 
A value of zero for \code{tol} means that the windows must be contiguous whereas negative values specify minimum overlaps.

If \code{sign!=NULL}, windows are only merged if they have the same sign of the log-FC and are not separated by intervening windows with opposite log-FC values. 
This can be useful to ensure consistent changes when summarizing adjacent DB regions.
However, it is not recommended for routine clustering in differential analyses as the resulting clusters will not be independent of the p-value.

Specification of \code{max.width} prevents the formation of excessively large clusters when many adjacent regions are present. 
Any cluster that is wider than \code{max.width} is split into multiple subclusters of (roughly) equal size.
Specifically, the cluster interval is partitioned into the smallest number of equally-sized subintervals where each subinterval is smaller than \code{max.width}.
Windows are then assigned to each subinterval based on the location of the window midpoints.
Suggested values range from 2000 to 10000 bp, but no limits are placed on the maximum size if it is \code{NULL}. 

The tolerance should reflect the minimum distance at which two regions of
enrichment are considered separate. If two windows are more than \code{tol}
apart, they \emph{will} be placed into separate clusters. In contrast, the
\code{max.width} value reflects the maximum distance at which two windows can be
considered part of the same region. 

Arbitrary regions can also be used in this function.
However, caution is required if any fully nested regions are present.
Clustering with \code{sign!=NULL} will throw an error as splitting by sign becomes undefined.
Splitting with \code{max.width!=NULL} will not fail, but cluster sizes may not be reduced if very large regions are present.

% For example, if you have a large region and lots of little nested regions within,
% the large region will be assigned to the middle subcluster whereas the nested elements
% will be assigned to their own subclusters. As a result, you end up with one giant 
% subcluster in which the other subclusters are nested.

If \code{ignore.strand=FALSE}, the entries in \code{regions} are split into their separate strands.
The function is run separately on the entries for each strand, and the results collated.
The \code{region} returned in the output will be stranded to reflect the strand of the contributing input regions.
This may be useful for strand-specific applications.

Note that, in the output, the cluster ID reported in \code{id} corresponds to the index of the cluster coordinates in the input \code{region}.
}

\examples{
x <- round(runif(10, 100, 1000))
gr <- GRanges(rep("chrA", 10), IRanges(x, x+40))
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
