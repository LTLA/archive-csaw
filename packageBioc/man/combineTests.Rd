\name{combineTests}
\alias{combineTests}

\title{Combine statistics across multiple tests}

\description{Combines p-values across clustered tests using Simes' method to control the cluster FDR.}

\usage{
combineTests(ids, tab, weight=rep(1, length(ids)))
}

\arguments{
	\item{ids}{an integer vector containing the cluster ID for each test}
	\item{tab}{a dataframe of results with \code{PValue}, \code{logCPM} and at least one \code{logFC} field for each test}
	\item{weight}{a numeric vector of weights for each window}
}

\value{
A dataframe with one row per cluster and the numeric fields \code{logCPM}, the
average log-CPM; \code{PValue}, the combined p-value; and \code{FDR}, the
q-value corresponding to the combined p-value. There is also one numeric field
representing the average of each supplied \code{logFC} in \code{tab}. The name
of each row corresponds to the sorted cluster IDs.
}

\details{
This function uses Simes' procedure to compute the combined p-value for each
cluster of tests. Each p-value represents evidence against the global null
hypothesis, i.e., all individual nulls are true in each cluster.  This may be
more relevant than examining each test individually when multiple tests in a
cluster represent parts of the same underlying event, e.g., genomic regions
consisting of clusters of windows.

Mean logFC and logCPM values are also computed across all tests in each
cluster. Multiple fields in \code{tab} containing the \code{logFC} substring
are allowed, e.g., to accommodate ANOVA-like contrasts. Note that the average
may not be a suitably informative metric when clusters are large and
heterogenous, in which case custom summaries may be required. The BH method is
also applied to control the FDR across all clusters.

The importance of each test within a cluster can be adjusted by supplying
different relative \code{weight} values. This may be useful for downweighting
low-confidence tests, e.g., those in repeat regions. In Simes' procedure,
weights are interpreted as relative frequencies of the tests in each cluster.
Note that these weights have no effect between clusters and will not be used to
adjust the computed FDR.

A simple clustering approach for windows is provided in
\code{\link{mergeWindows}}. However, anything can be used so long as it does
not compromise type I error control, e.g., promoters, gene bodies,
independently called peaks. 
} 

\examples{ 
ids <- round(runif(100, 1, 10))
tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100),
	    PValue=rbeta(100, 1, 2))
combined <- combineTests(ids, tab)
head(combined)

# With window weighting.
w <- round(runif(100, 1, 5))
combined <- combineTests(ids, tab, weight=w)
head(combined)

# With multiple log-FCs.
tab$logFC.whee <- rnorm(100, 5)
combined <- combineTests(ids, tab)
head(combined)
}

\seealso{
\code{\link{mergeWindows}}
}

\author{Aaron Lun}

\references{
Simes RJ (1986). An improved Bonferroni procedure for multiple tests of
significance. \emph{Biometrika} 73, 751-754.

Benjamini Y and Hochberg Y (1995). Controlling the false discovery rate: a
practical and powerful approach to multiple testing. \emph{J. R. Stat. Soc.
Series B} 57, 289-300. 

Lun ATL and Smyth GK (2014). De novo detection of differentially bound regions
for ChIP-seq data using peaks and windows: controlling error rates correctly.
\emph{Nucleic Acids Res.} 42, e95
}

\keyword{testing}
