\name{empiricalFDR}
\alias{empiricalFDR}

\title{Control the empirical FDR}

\description{Control the empirical FDR across clusters for comparisons to negative controls, based on tests that are significant in the wrong direction.}

\usage{
empiricalFDR(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL, neg.down=TRUE) 
}

\arguments{
\item{ids}{an integer vector containing the cluster ID for each test}
\item{tab}{a dataframe of results with \code{PValue} and at least one \code{logFC} field for each test}
\item{weight}{a numeric vector of weights for each window, defaults to 1 for each test}
\item{pval.col}{an integer scalar or string specifying the column of \code{tab} containing the p-values}
\item{fc.col}{an integer scalar or string specifying the column of \code{tab} containing the log-fold changes}
\item{neg.down}{a logical scalar indicating if negative log-fold changes correspond to the \dQuote{wrong} direction}
}

\details{
Some experiments involve comparisons to negative controls where there should be no signal/binding.
In such case, genuine differences should only occur in one direction, i.e., up in the non-control samples.
Thus, the number of significant tests that change in the wrong direction can be used as an estimate of the number of false positives.

This function converts two-sided p-values in \code{tab[,pval.col]} to one-sided counterparts in the wrong direction.
It combines the one-sided p-values for each cluster using \code{\link{combineTests}}.
The number of significant clusters at some p-value threshold represents the estimated number of false positive clusters.

The same approach is applied for one-sided p-values in the right direction, where the number of detected clusters at the threshold represents the total number of discoveries.
Dividing the number of false positives by the nmber of discoveries yields the empirical FDR at each p-value threshold.
Monotonicity is enforced (i.e., the empirical FDR only decreases with decreasing p-value) as is the fact that the empirical FDR must be below unity.

The p-values specified in \code{pval.col} are assumed to be originally computed from some two-sided test,
where the direction of change is independent of the magnitude of the p-value under the null hypothesis.
This rules out p-values computed from, e.g., ANODEV where multiple contrasts are tested at once.

Control of the empirical FDR is best used for very noisy data sets where the BH method is not adequate.
While the BH method protects against statistical false positives, the empirical FDR also protects against experimental false positives, e.g., due to non-specific binding.
The downside is that the empirical FDR calculation relies on the availability of a good estimate of the number of false positives.
The BH method in \code{\link{combineTests}} is more statistically rigorous and should be preferred for routine analyses.
}

\value{
A data frame containing one row per cluster, with various fields:
\itemize{
\item A numeric field containing the one-sided p-value for each cluster in the right direction.
This field is named \code{PValue} if \code{pval.col=NULL}, otherwise its name is set to \code{colnames(tab[,pval.col]}.
\item A numeric field \code{FDR}, containing the empirical FDR corresponding to the p-value threshold equal to the value in \code{PValue}.
}
All other fields are the same as those returned by \code{\link{combineTests}}.
The exception is the \code{direction} field, which is not returned as it is not informative for one-sided tests.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{combineTests}}
}

\references{
Zhang Y, Liu T, Meyer CA et al. (2008). Model-based Analysis of ChIP-Seq (MACS). \emph{Genome Biol.} 9, R137.
}

\examples{
ids <- round(runif(100, 1, 10))
tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
empirical <- empiricalFDR(ids, tab)
head(empirical)
}
