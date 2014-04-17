\name{getBestTest}
\alias{getBestTest}

\title{Get the best test in a cluster}

\description{Find the test with the strongest evidence for rejection of the null in each cluster}

\usage{
getBestTest(ids, tab, weight=rep(1, length(ids)))
}

\arguments{
	\item{ids}{an integer vector containing the cluster ID for each test}
	\item{tab}{a table of results with a \code{PValue} field for each test}
	\item{weight}{a numeric vector of weights for each window}
}

\value{
A dataframe with one row per cluster and the numeric fields \code{best}, the index for the best
test in the cluster; \code{PValue}, the adjusted p-value for that test; and \code{FDR}, the 
q-value corresponding to the adjusted p-value.
}

\details{
This function identifies the test with the lowest p-value as that with the strongest evidence against 
the null in each cluster. The p-value of the chosen test is adjusted using the Bonferroni correction,
based on the total number of tests in the parent cluster. The procedure is necessary to obtain
strong control of the FWER across all configurations of the nulls within each cluster. In particular,
this means that the best test can be taken from each cluster for further consideration.

% i.e. The configuration, in this case, is taking the best test.

The importance of each window in each cluster can be adjusted by supplying different relative 
\code{weight} values. Each weight is interpreted as a different threshold for each test in the cluster.
Larger weights correspond to lower thresholds, i.e., less evidence is needed to reject the null for
tests deemed to be more important. This may be useful for upweighting particular tests, e.g., windows
containing a motif for the TF of interest.
}

\examples{ 
ids<-round(runif(100, 1, 10))
tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100),
	PValue=rbeta(100, 1, 2))
combined<-getBestTest(ids, tab)
head(combined)

# With window weighting.
w<-round(runif(100, 1, 5))
combined<-getBestTest(ids, tab, weight=w)
head(combined)
}

\seealso{
\code{\link{combineTests}}, \code{\link{mergeWindows}}, 
}

\author{Aaron Lun}

\references{
Genovese CR, Roeder K and Wasserman L (2006). False discovery control with p-value weighting. \emph{Biometrika} 93, 509-524.
}
