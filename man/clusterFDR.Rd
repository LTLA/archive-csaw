\name{clusterFDR}
\alias{clusterFDR}
\alias{controlClusterFDR}

\title{Compute the cluster-level FDR}
\description{Compute the FDR across clusters based on the test-level FDR threshold}

\usage{
clusterFDR(ids, threshold, weight=NULL)
controlClusterFDR(target, adjp, FUN, ..., weight=NULL, grid.param=NULL)
}

\arguments{
\item{ids}{an integer vector of cluster IDs for each significant test below \code{threshold}}
\item{threshold}{a numeric scalar, specifying the FDR threshold used to define the significant tests}
\item{target}{a numeric scalar specifying the desired cluster-level FDR threshold}
\item{adjp}{a numeric vector of window-level adjusted p-values}
\item{FUN}{a clustering function that takes a logical vector indicating which windows are significant, and returns an integer vector of cluster IDs (see below)}
\item{...}{additional arguments to be passed to \code{FUN}}
\item{weight}{a numeric vector of frequency weights, for internal use}
\item{grid.param}{a named list of grid search parameters, see Details}
} 

\section{Definition of the cluster-level FDR}{
The \code{clusterFDR} function computes an informal estimate of the cluster-level FDR, where each cluster is formed by aggregating only significant tests.
In the context of ChIP-seq, each significant test refers to a DB window that is detected at a FDR below \code{threshold}.
The idea is to obtain an error rate while reporting the precise coordinates of a DB subinterval in a complex region.

The cluster-level FDR is defined as the proportion of reported clusters that have no true positives.
Simply using \code{threshold} is not appropriate, as the cluster- and window-level FDRs are not equivalent. 
This function also differs from the standard pipeline that is based on \code{\link{combineTests}}.
Specifically, region definition in \code{\link{combineTests}} must be independent of DB so that precise coordinates of the DB subinterval cannot be reported.
This is overcome here, by clustering directly on DB windows and applying post-hoc control of the cluster-level FDR.

Note that the calculation of the cluster-level FDR here is not statistically rigorous.
In particular, the observed number of false positive tests is estimated based on \code{threshold} and the total number of significant tests.
This is not guaranteed to be an upper bound, especially with few or correlated tests. 
Thus, users should use the standard \code{\link{combineTests}}-based pipeline wherever possible.
Clustering on significant windows should only be performed where the precise coordinates of the DB subinterval are important for interpretation.
}

\section{Searching for the best threshold}{
\code{controlClusterFDR} will identify the window-level FDR threshold required to control the cluster-level FDR at \code{target}.
The former is not a simple function of the latter (neither continuous nor guaranteed to be monotonic), so a grid search is used.
Clusters of significant windows are identified at each window-level threshold, and the corresponding cluster-level FDR is computed with \code{clusterFDR}.

At each iteration, the grid point with the closest cluster-level FDR to \code{target} is chosen and the grid is recentered around that point.
The size of the grid is also scaled down to provide greater resolution.
Users can tune the settings of the grid search by specifying elements in \code{grid.param} as:
\describe{
\item{\code{length}:}{an integer scalar specifying the length of the grid, defaults to 21}
\item{\code{range}:}{a numeric scalar indicating the range of the grid in logit units, defaults to 20}
\item{\code{iter}:}{an integer scalar indicating the number of iterations of the grid search}
\item{\code{scale}:}{a numeric scalar specifying how the range of the grid should be downscaled at each iteration}
}

The \code{FUN} argument should be a function that accepts a logical vector specifying significance, and returns an integer vector of cluster IDs.
If, for example, it accepts an input vector \code{ix}, then the output should contain cluster IDs corresponding to the entries of \code{which(ix)}.
This is because cluster IDs are only defined for significant tests, given that only those tests are used for clustering.

An additional requirement is that the returned window-level FDR threshold should be less than \code{target}.
Thus, each window should be significantly DB on its own merits before it is placed into a cluster.
This protects against scenarios where very large thresholds yield low cluster-level FDRs, due to the formation of a few large clusters.
}

\section{Note about weights}{ 
In both functions, the \code{weight} argument is assumed to contain frequency weights of significant tests/windows.
For example, a weight of 2 for a test would be equivalent to repeating that test (i.e., repeating the same window so it shows up twice in your analysis).
These weights should be the same as those used during weighted FDR control to compute adjusted p-values.
In general, you should not set this argument unless you know what you're doing.
}

\value{
For \code{clusterFDR}, a numeric scalar is returned as the cluster-level FDR.

For \code{controlClusterFDR}, a list is returned containing two numeric scalars -- \code{threshold}, the window-level FDR threshold to control the cluster-level FDR near \code{target}; and \code{FDR}, the estimate of the cluster-level FDR corresponding to \code{threshold}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{mergeWindows}},
\code{\link{combineTests}},
\code{\link{clusterWindows}}
}

\examples{
# Setting up the windows and p-values.
set.seed(100)
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
