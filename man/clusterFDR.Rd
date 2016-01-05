\name{clusterFDR}
\alias{clusterFDR}
\alias{controlClusterFDR}

\title{Compute the cluster-level FDR}
\description{Compute the FDR across clusters based on the test-level FDR threshold}

\usage{
clusterFDR(ids, threshold)
controlClusterFDR(target, adjp, FUN, ..., maxiter=5, contract=4)
}

\arguments{
\item{ids}{an integer vector of cluster IDs for each significant test below \code{threshold}}
\item{threshold}{a numeric scalar, specifying the FDR threshold used to define the significant tests}
\item{target}{a numeric scalar specifying the desired cluster-level FDR threshold}
\item{adjp}{a numeric vector of window-level adjusted p-values}
\item{FUN}{a clustering function that takes a logical vector indicating which windows are significant, and returns an integer vector of cluster IDs}
\item{...}{additional arguments to be passed to \code{FUN}}
\item{maxiter}{an integer scalar indicating the maximum number of iterations for the grid search}
\item{contract}{a numeric scalar specifying the fold decrease in the grid range at every iteration}
} 

\details{
The \code{clusterFDR} function computes an informal estimate of the cluster-level FDR, where each cluster is formed by aggregating only significant tests.
In the context of ChIP-seq, each significant test refers to a DB window that is detected at a FDR below \code{threshold}.
The idea is to obtain an error rate while reporting the precise coordinates of a DB subinterval in a complex region.

The cluster-level FDR is defined as the proportion of reported clusters that have no true positives.
Simply using \code{threshold} is not appropriate, as the cluster- and window-level FDRs are not equivalent. 
This function also differs from the standard pipeline that is based on \code{\link{combineTests}}.
Specifically, region definition in the standard pipeline must be independent of DB.
Precise coordinates of the DB subinterval cannot be reported.

Users should note that the calculation of the cluster-level FDR here is not statistically rigorous.
In particular, the observed number of false positive tests is estimated based on \code{threshold} and the total number of significant tests.
This is not guaranteed to be an upper bound, especially when the observed window-level FDR is variable.
Thus, users should use the standard \code{\link{combineTests}}-based pipeline wherever possible.
Clustering on significant windows should only be performed where the precise coordinates of the DB subinterval are important for interpretation.

The \code{controlClusterFDR} will identify the window-level FDR threshold required to control the cluster-level FDR at \code{target}.
This is necessary as the latter is not a simple function of the former.
A grid search is used to identify clusters of significant windows at each window-level threshold, and the corresponding cluster-level FDR is computed with \code{clusterFDR}.
Users should note that the discontinuous nature of the cluster-level FDR means that the returned value may not be below \code{target}.
}

\value{
For \code{clusterFDR}, a numeric scalar is returned as the cluster-level FDR.

For \code{controlClusterFDR}, a list is returned containing two numeric scalars -- \code{threshold}, the window-level FDR threshold to control the cluster-level FDR near \code{target}; and \code{cluster.FDR}, the estimate of the cluster-level FDR corresponding to \code{threshold}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{mergeWindows}},
\code{\link{combineTests}}
}

\examples{
# Setting up the windows and p-values.
windows <- GRanges("chrA", IRanges(1:1000, 1:1000))
test.p <- runif(1000)
test.p[c(1:10, 100:110, 220:240)] <- 0 # 3 significant subintervals.

# Defining significant windows.
threshold <- 0.05
is.sig <- p.adjust(test.p, method="BH") <= threshold

# Assuming that we only cluster significant windows.
merged <- mergeWindows(windows[is.sig], tol=0)
clusterFDR(merged$id, threshold)

# Setting up another example with more subintervals.
test.p <- runif(1000)
test.p[rep(1:2, 50) + rep(0:49, each=2) * 20] <- 0 
adj.p <- p.adjust(test.p, method="BH")
clusterFUN <- function(x) { mergeWindows(windows[x], tol=0)$id }
controlClusterFDR(0.05, adj.p, clusterFUN)             
}

\keyword{testing}
