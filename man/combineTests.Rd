\name{combineTests}
\alias{combineTests}

\title{Combine statistics across multiple tests}

\description{Combines p-values across clustered tests using Simes' method to control the cluster FDR.}

\usage{
combineTests(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL)
}

\arguments{
	\item{ids}{an integer vector or factor containing the cluster ID for each test}
	\item{tab}{a dataframe of results with \code{PValue} and at least one \code{logFC} field for each test}
	\item{weight}{a numeric vector of weights for each window, defaults to 1 for each test}
	\item{pval.col}{an integer scalar or string specifying the column of \code{tab} containing the p-values}
	\item{fc.col}{an integer or character vector specifying the columns of \code{tab} containing the log-fold changes}
}

\value{
A dataframe with one row per cluster and various fields:
\itemize{
\item An integer field \code{nWindows}, specifying the total number of windows in each cluster.
\item Two integer fields \code{*.up} and \code{*.down} for each log-FC column in \code{tab}, containing the number of windows with log-FCs above 0.5 or below -0.5, respectively.
\item A numeric field containing the combined p-value. 
If \code{pval.col=NULL}, this column is named \code{PValue}, otherwise its name is set to \code{colnames(tab[,pval.col])}.
\item A numeric field \code{FDR}, containing the q-value corresponding to the combined p-value.
\item A character field \code{direction} (if \code{fc.col} is of length 1), specifying the dominant direction of change for windows in each cluster.
}
Each row is named according to the ID of the corresponding cluster.
}

\details{
This function uses Simes' procedure to compute the combined p-value for each cluster of tests with the same value of \code{ids}.
Each combined p-value represents evidence against the global null hypothesis, i.e., all individual nulls are true in each cluster. 
This may be more relevant than examining each test individually when multiple tests in a cluster represent parts of the same underlying event, e.g., genomic regions consisting of clusters of windows.
The BH method is also applied to control the FDR across all clusters.

The importance of each test within a cluster can be adjusted by supplying different relative \code{weight} values. 
This may be useful for downweighting low-confidence tests, e.g., those in repeat regions. 
In Simes' procedure, weights are interpreted as relative frequencies of the tests in each cluster. 
Note that these weights have no effect between clusters and will not be used to adjust the computed FDR.

By default, the relevant fields in \code{tab} are identified by matching the column names to their expected values.
Multiple fields in \code{tab} containing the \code{logFC} substring are allowed, e.g., to accommodate ANOVA-like contrasts. 
The p-value column is expected to be named as \code{PValue}.
If the column names are different from what is expected, specification of the correct columns can be performed using \code{pval.col} and \code{fc.col}.
This will overwrite any internal selection of the appropriate fields.

A simple clustering approach for windows is provided in \code{\link{mergeWindows}}. 
However, anything can be used so long as it is independent of the p-values and does not compromise type I error control, e.g., promoters, gene bodies, independently called peaks. 
Any tests with \code{NA} values for \code{ids} will be ignored.
} 

\section{Determining the direction of DB}{
This function will report the number of windows with log-fold changes above 0.5 and below -0.5, to give some indication of whether binding increases or decreases in the cluster.
If a cluster contains non-negligble numbers of \code{up} and \code{down} windows, this indicates that there may be a complex DB event within that cluster.
Similarly, complex DB may be present if the total number of windows is larger than the number of windows in either category (i.e., change is not consistent across the cluster).
Note that the threshold of 0.5 is arbitrary and has no impact on the significance calculations.

When only one log-fold change column is specified, \code{combineTests} will determine which DB direction contributes to the combined p-value.
This is done by considering whether the combined p-value would increase if all tests in one direction were assigned p-values of unity.
If there is an increase, then tests changing in that direction must contribute to the calculations in Simes' method.
In this manner, clusters are labelled based on whether they are driven by tests with positive (\code{"up"}) or negative log-fold changes (\code{"down"}) or both (\code{"mixed"}).

The label for each cluster is stored as the \code{direction} field in the returned data frame.
However, keep in mind that the label only describes the direction of change among the most significant tests in the cluster.
Clusters with complex DB may still be labelled as changing in only one direction, 
if the tests changing in one direction have much lower p-values than the tests changing in the other direction (even if both sets of p-values are significant).
}

\examples{ 
ids <- round(runif(100, 1, 10))
tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
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

# Manual specification of column IDs.
combined <- combineTests(ids, tab, fc.col=c(1,4), pval.col=3)
head(combined)

combined <- combineTests(ids, tab, fc.col="logFC.whee", pval.col="PValue")
head(combined)
}

\seealso{
\code{\link{mergeWindows}}
}

\author{Aaron Lun}

\references{
Simes RJ (1986). An improved Bonferroni procedure for multiple tests of significance. \emph{Biometrika} 73, 751-754.

Benjamini Y and Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{J. R. Stat. Soc. Series B} 57, 289-300. 

Benjamini Y and Hochberg Y (1997). Multiple hypotheses testing with weights. \emph{Scand. J. Stat.} 24, 407-418.

Lun ATL and Smyth GK (2014). De novo detection of differentially bound regions for ChIP-seq data using peaks and windows: controlling error rates correctly. \emph{Nucleic Acids Res.} 42, e95
}

\keyword{testing}
