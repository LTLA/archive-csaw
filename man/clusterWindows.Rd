\name{clusterWindows}
\alias{clusterWindows}

\title{Cluster DB windows into clusters}

\description{Clusters significant windows into clusters while controlling the cluster-level FDR.}

\usage{
clusterWindows(regions, tab, target=0.05, pval.col=NULL, 
    tol, sign=NULL, ..., grid.param=NULL) 
}

\arguments{
	\item{regions}{a GRanges or RangedSummarizedExperiment object containing window coordinates}
    \item{tab}{a dataframe of results with a \code{PValue} field for each window}
    \item{target}{a numeric scalar indicating the desired cluster-level FDR}
    \item{pval.col}{an argument specifying the column of \code{tab} with the p-values, see \code{\link{combineTests}} for details}
	\item{tol, sign, ...}{arguments to be passed to \code{\link{mergeWindows}}}
    \item{grid.param}{arguments to be passed to \code{\link{controlClusterFDR}}}
}

\value{
A named list similar to that reported by \code{\link{mergeWindows}} with an ID vector in \code{id} and region coordinates of each cluster in \code{region}.
Non-significant windows are marked with \code{NA} values in \code{ids}.
An additional element \code{FDR} is also included, representing the estimate of the cluster-level FDR for the returned regions.
}

\details{
Windows are identified as DB based on the adjusted p-values in \code{tab}.
Only these DB windows are then used directly for clustering via \code{\link{mergeWindows}}, in order to identify DB regions.
This is not blind to the DB status so post-hoc control of the cluster-level FDR is applied by using \code{\link{controlClusterFDR}}.
The aim is to provide some interpretable results when DB-blind clustering is not appropriate, e.g., for diffuse marks involving long stretches of the genome.
Reporting each marked stretch in its entirety would be cumbersome, so this method allows the DB subintervals to be identified directly.
}

\examples{
set.seed(10)
x <- round(runif(100, 100, 1000))
gr <- GRanges("chrA", IRanges(x, x+5))
tab <- data.frame(PValue=rbeta(length(x), 1, 50))

clusterWindows(gr, tab, tol=10)
clusterWindows(gr, tab, tol=10, sign=rbinom(100, 1, 0.5)==0L)
}

\seealso{
\code{\link{mergeWindows}},
\code{\link{controlClusterFDR}}
}

\author{Aaron Lun}

\keyword{clustering}
