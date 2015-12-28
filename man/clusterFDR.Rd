\name{clusterFDR}
\alias{clusterFDR}

\title{Compute the cluster-level FDR}
\description{Compute the FDR across clusters based on the test-level FDR threshold}

\usage{
clusterFDR(ids, threshold)
}

\arguments{
\item{ids}{an integer vector of cluster IDs for each significant test below \code{threshold}}
\item{threshold}{a numeric scalar, specifying the FDR threshold used to define the significant tests}
} 

\details{
This function computes an informal estimate of the cluster-level FDR, where each cluster is formed by aggregating only significant tests.
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

In conclusion, users should use the standard \code{\link{combineTests}}-based pipeline wherever possible.
Clustering on significant windows should only be performed where the precise coordinates of the DB subinterval are important for interpretation.
}

\value{
A numeric scalar as the cluster-level FDR.
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
}

\keyword{testing}
