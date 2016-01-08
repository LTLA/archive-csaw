\name{consolidateClusters}
\alias{consolidateClusters}

\title{Consolidate DB clusters}
\description{Consolidate DB results from multiple analyses with cluster-level FDR control.}

\usage{
consolidateClusters(data.list, result.list, equiweight=TRUE, ...)
}

\arguments{
\item{data.list}{a list of RangedSummarizedExperiment and/or GRanges objects}
\item{result.list}{a list of data frames containing the DB test results for each entry of \code{data.list}}
\item{equiweight}{a logical scalar indicating whether equal weighting from each analysis should be enforced}
\item{...}{arguments to be passed to \code{\link{clusterWindows}}}
}

\value{
A named list is returned containing \code{id}, a list of integer vectors indicating the cluster ID for each window in \code{data.list}; \code{region}, a GRanges object containing the coordinates for each cluster; and \code{FDR}, a numeric scalar containing the cluster-level FDR estimate.
}

\details{
This function consolidates DB results from multiple analyses, typically involving different window sizes.
The aim is to provide comprehensive detection of DB at a range of spatial resolutions.
Significant windows from each analysis are identified and used for clustering with \code{\link{clusterWindows}}.
This represents the post-hoc counterpart to \code{\link{consolidateSizes}}.

Some effort is required to equalize the contribution of the results from each analysis.
This is done by setting \code{equiweight=TRUE}, where the weight of each window is inversely proportional to the number of windows from that analysis.
These weights are used as frequency weights for window-level FDR control (to identify DB windows prior to clustering).
Otherwise, the final results would be dominated by large number of small windows.

% Fully nested windows of opposite sign will break it, so don't try it.
}

\seealso{
\code{\link{clusterWindows}}
}

\author{
Aaron Lun
}

\examples{
# Making up some GRanges.
low <- GRanges("chrA", IRanges(runif(100, 1, 1000), width=5))
med <- GRanges("chrA", IRanges(runif(40, 1, 1000), width=10))
high <- GRanges("chrA", IRanges(runif(10, 1, 1000), width=20))

# Making up some DB results.
dbl <- data.frame(logFC=rnorm(length(low)), PValue=rbeta(length(low), 1, 20))
dbm <- data.frame(logFC=rnorm(length(med)), PValue=rbeta(length(med), 1, 20))
dbh <- data.frame(logFC=rnorm(length(high)), PValue=rbeta(length(high), 1, 20))
result.list <- list(dbl, dbm, dbh)

# Consolidating.
cons <- consolidateClusters(list(low, med, high), result.list, tol=20)
cons$region
cons$id
cons$FDR

# Without weights.
cons <- consolidateClusters(list(low, med, high), result.list,
            tol=10, equiweight=FALSE)
cons$region
cons$id
cons$FDR
}

\keyword{clustering}
