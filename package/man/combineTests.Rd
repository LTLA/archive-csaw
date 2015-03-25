\name{combineTests}
\alias{combineTests}

\title{Combine statistics across multiple tests}

\description{Combines p-values across clustered tests using Simes' method to control the cluster FDR.}

\usage{
combineTests(ids, tab, weight=NULL, pval.col=NULL, fc.col=NULL)
}

\arguments{
	\item{ids}{an integer vector containing the cluster ID for each test}
	\item{tab}{a dataframe of results with \code{PValue}, \code{logCPM} and at least one \code{logFC} field for each test}
	\item{weight}{a numeric vector of weights for each window, defaults to 1 for each test}
	\item{pval.col}{an integer scalar specifying the column of \code{tab} containing the p-values, or a character string containing the name of that column}
	\item{fc.col}{an integer vector specifying the columns of \code{tab} containing the log-fold changes, or a character vector containing the names of those columns}
}

\value{
A dataframe with one row per cluster and the numeric fields \code{PValue}, the combined p-value; and \code{FDR}, the q-value corresponding to the combined p-value.
An integer field \code{nWindows} specifies the total number of windows in each cluster.
There are also two integer fields \code{*.up} and \code{*.down} for each log-FC column in \code{tab}, containing the number of windows with log-FCs above 0.5 or below -0.5, respectively.
The name of each row represents the ID of the corresponding cluster.
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
If the column names are different from what is expected, specification of the correct columns can be performed using \code{pval.col} and \code{fc.col}.
This will overwrite any internal selection of the appropriate fields.

This function will report the number of windows with log-fold changes above 0.5 and below -0.5, to give some indication of whether binding increases or decreases in the cluster.
If a cluster contains non-negligble numbers of \code{up} and \code{down} windows, this indicates that there may be a complex DB event within that cluster.
Similarly, complex DB may be present if the total number of windows is larger than the number of windows in either category (i.e., change is not consistent across the cluster).
Note that the threshold of 0.5 is arbitrary and has no impact on the significance calculations.

A simple clustering approach for windows is provided in \code{\link{mergeWindows}}. 
However, anything can be used so long as it does not compromise type I error control, e.g., promoters, gene bodies, independently called peaks. 
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
