\name{mixedClusters}
\alias{mixedClusters}

\title{Tests for mixed DB clusters}

\description{Intersects two one-sided tests to determine if a cluster contains DB events in both directions.}

\usage{
mixedClusters(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL) 
}

\arguments{
\item{ids}{an integer vector or factor containing the cluster ID for each test}
\item{tab}{a dataframe of results with \code{PValue} and at least one \code{logFC} field for each test}
\item{weight}{a numeric vector of weights for each window, defaults to 1 for each test}
\item{pval.col}{an integer scalar or string specifying the column of \code{tab} containing the p-values}
\item{fc.col}{an integer scalar or string specifying the columns of \code{tab} containing the log-fold changes}
}

\value{
A data frame containing one row per cluster, with various fields:
\itemize{
\item A numeric field containing the IUT p-value for each cluster.
This field is named \code{PValue} if \code{pval.col=NULL}, otherwise its name is set to \code{colnames(tab[,pval.col]}.
\item A numeric field \code{FDR}, containing the adjusted p-value for each cluster.
}
All other fields are the same as those returned by \code{\link{combineTests}}.
The exception is the \code{direction} field, which is not returned as the test is done explicitly for mixed clusters.
}

\details{
This function converts two-sided p-values to one-sided counterparts for each direction of log-fold change.
For each direction, the corresponding one-sided p-values are combined to yield a combined p-value for each cluster.
Each cluster is associated with two combined p-values (one in each direction), which are intersected using the Berger's intersection-union test (IUT).
This tests against the null hypothesis that either direction is not significant.
In short, a low p-value is only possible if there is DB in both directions.
This formally identifies mixed clusters corresponding to complex DB events.
}

\seealso{
\code{\link{combineTests}}
}

\author{Aaron Lun}

\examples{
ids <- round(runif(100, 1, 10))
tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
mixed <- mixedClusters(ids, tab)
head(mixed)
}
