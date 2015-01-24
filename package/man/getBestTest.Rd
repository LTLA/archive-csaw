\name{getBestTest}
\alias{getBestTest}

\title{Get the best test in a cluster}

\description{Find the test with the strongest evidence for rejection of the null in each cluster.}

\usage{
getBestTest(ids, tab, by.pval=TRUE, weight=rep(1, length(ids)), pval.col=NULL, cpm.col=NULL)
}

\arguments{
	\item{ids}{an integer vector containing the cluster ID for each test}
	\item{tab}{a table of results with a \code{PValue} field for each test}
	\item{by.pval}{a logical scalar, indicating whether selection should be performed on corrected p-values}
	\item{weight}{a numeric vector of weights for each window}
	\item{pval.col}{an integer scalar specifying the column of \code{tab} containing the p-values}
	\item{cpm.col}{an integer scalar specifying the column of \code{tab} containing the log-CPM values}
}

\value{
A dataframe with one row per cluster and the numeric fields \code{best}, the index for the best test in the cluster; \code{PValue}, the adjusted p-value for that test; and \code{FDR}, the q-value corresponding to the adjusted p-value.
Note that the p-value column may be named differently if \code{pval.col} is specified.
Cluster IDs are stored as the row names.
}

\details{
If \code{by.pval=TRUE}, this function identifies the test with the lowest p-value as that with the strongest evidence against the null in each cluster.
The p-value of the chosen test is adjusted using the Bonferroni correction, based on the total number of tests in the parent cluster. 
This is necessary to obtain strong control of the family-wise error rate such that the best test can be taken from each cluster for further consideration.

% i.e. The configuration, in this case, is taking the best test.

The importance of each window in each cluster can be adjusted by supplying different relative \code{weight} values. 
Each weight is interpreted as a different threshold for each test in the cluster.
Larger weights correspond to lower thresholds, i.e., less evidence is needed to reject the null for tests deemed to be more important. 
This may be useful for upweighting particular tests, e.g., windows containing a motif for the TF of interest.

Note the difference between this function and \code{\link{combineTests}}. 
The latter presents evidence for any rejections within a cluster. 
This function specifies the exact location of the rejection in the cluster, which may be more useful in some cases but at the cost of conservativeness. 
In both cases, clustering procedures such as \code{\link{mergeWindows}} can be used to identify the cluster.

% The vagueness of combineTests may be good enough in most applications (i.e. wanting to get a location
% to look at the genomic context, or in instances where differential binding is obvious). If error control
% at specific locations is needed, then getBestTests is probably more appropriate..

If \code{by.pval=FALSE}, the best test is defined as that with the highest log-CPM value. 
This should be independent of the p-value so no adjustment is necessary. Weights are not applied here. 
This mode may be useful when abundance is correlated to rejection under the alternative hypothesis, e.g., picking high-abundance regions that are more likely to contain peaks.

By default, the relevant fields in \code{tab} are identified by matching the column names to their expected values.
If the column names are different from what is expected, specification of the correct column can be performed using \code{pval.col} and \code{cpm.col}
}

\examples{ 
ids <- round(runif(100, 1, 10))
tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), 
    PValue=rbeta(100, 1, 2))
combined <- getBestTest(ids, tab)
head(combined)

# With window weighting.
w <- round(runif(100, 1, 5))
combined <- getBestTest(ids, tab, weight=w)
head(combined)

# By logCPM.
combined <- getBestTest(ids, tab, by.pval=FALSE)
head(combined)
}

\seealso{
\code{\link{combineTests}}, \code{\link{mergeWindows}} 
}

\author{Aaron Lun}

\references{
Genovese CR, Roeder K and Wasserman L (2006). False discovery control with
p-value weighting. \emph{Biometrika} 93, 509-524.
}

\keyword{testing}
