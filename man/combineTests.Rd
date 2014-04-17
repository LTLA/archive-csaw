\name{combineTests}
\alias{combineTests}

\title{Combine statistics across multiple tests}

\description{Combines p-values across clustered tests using Simes' method to control the cluster FDR.}

\usage{
combineTests(ids, tab, weight=rep(1, length(ids)))
}

\arguments{
	\item{ids}{an integer vector containing the cluster ID for each test}
	\item{tab}{a table of results with \code{PValue}, \code{logCPM} and at least one \code{logFC} field for each test}
	\item{weight}{a numeric vector of weights for each window}
}

\value{
A dataframe with one row per cluster and the numeric fields \code{logFC}, the average of each supplied log-FC
across the cluster; \code{logCPM}, the average log-CPM; \code{PValue}, the combined p-value; 
and \code{FDR}, the q-value corresponding to the combined p-value. The name of each row corresponds
to the sorted cluster IDs.
}

\details{
This function uses Simes' procedure to compute the combined p-value for each cluster of tests. Each p-value 
represents evidence against the global null hypothesis, i.e., all individual nulls are true in each cluster.
This may be more relevant than examining each test individually when multiple tests in a cluster represent 
parts of the same underlying event, e.g., genomic regions consisting of clusters of windows.

Mean logFC and logCPM values are also computed across all tests in each cluster. Multiple fields in \code{tab} 
containing the \code{logFC} pattern are allowed, e.g., to accommodate ANOVA-like contrasts. Note that the 
average may not be a suitably informative metric when clusters are large and heterogenous, so custom summaries 
may be required. The BH method is also applied to control the FDR across all clusters.

The importance of each test within a cluster can be adjusted by supplying different relative \code{weight} 
values. This may be useful for downweighting low-confidence tests, e.g., those in repeat regions. In Simes' 
procedure, weights are interpreted as relative frequencies of the tests in each cluster. Note that these weights 
have no effect between clusters and will not be used to adjust the computed FDR.

A simple clustering approach for windows is provided in \code{\link{mergeWindows}}. However, anything can be 
used so long as it does not compromise type I error control, e.g., promoters, gene bodies, independently called 
peaks. In all cases, consecutive cluster IDs (i.e., within \code{[1, N]} for \code{N} clusters) are recommended. 
This allows the direct use of \code{ids} as an index for the output table.
} 

\examples{ 
counts<-matrix(rnbinom(400, mu=10, size=20), ncol=4)
y<-DGEList(counts)
design<-model.matrix(~factor(c("a", "a", "b", "b")))
fit<-glmFit(y, design, dispersion=0.05)
out<-glmLRT(fit)
ids<-round(runif(100, 1, 10))
combined<-combineTests(ids, out$table)
head(combined)

# With window weighting.
w<-round(runif(100, 1, 5))
combined<-combineTests(ids, out$table, weight=w)
head(combined)

# With cluster weighting, defined as the mean of window weights.
cluster.w<-as.numeric(by(w, list(ids), FUN=mean))
o<-order(combined$PValue)
total.w<-sum(cluster.w)
cs.w<-cumsum(cluster.w[o])
combined$FDR[o]<-rev(cummin(rev(combined$PValue[o]*total.w/cs.w)))
}

\seealso{
\code{\link{mergeWindows}}
}

\author{Aaron Lun}

\references{
Simes RJ (1986). An improved Bonferroni procedure for multiple tests of significance. \emph{Biometrika} 73, 751-754.

Benjamini Y and Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society Series B} 57, 289-300. 
}
