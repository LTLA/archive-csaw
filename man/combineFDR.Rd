\name{combineFDR}
\alias{combineFDR}

\title{Combine statistics across multiple windows}

\description{Combines p-values across clustered windows using Simes' method to control the cluster FDR.}

\usage{
combineFDR(ids, tab, weight=rep(1, length(ids)))
}

\arguments{
	\item{ids}{an integer vector containing the cluster ID for each window}
	\item{tab}{a table of results with p-values, log-FC and log-CPM values for each window}
	\item{weight}{a numeric vector of weights for each window}
}

\value{
A dataframe with one row per cluster and the numeric fields \code{logFC}, the average log-FC across
all windows in each cluster; \code{logCPM}, the average log-CPM; \code{PValue}, the combined p-value; 
and \code{FDR}, the q-value corresponding to the combined p-value.
}

\details{
Reporting the window-based FDR is often misleading as interpretation and validation of ChIP-seq binding
sites is performed on a region-by-region basis. This means that applying the BH method across p-values for
all windows will be optimistic unless one plans to validate every overlapping DE window. 

In a window-based framework, the more relevant FDR is that across clusters of windows. This can be 
controlled by computing the p-value for the global null hypothesis across a cluster i.e. no enrichment
at any window within the cluster. Simes' method can be used to compute a combined p-value for this purpose.

The cluster FDR can then be simply controlled by applying the BH method to the combined p-values for all clusters.
Rejection of the global null indicates that the cluster is worth further investigation. The use of window-based counts to 
compute the original p-values provides substantial robustness to changes in the definition of a cluster.

A simple clustering approach is provided in \code{mergeWindows}. However, anything can be used so long as it 
does not compromise type I error control. For example, peak-calling programs can be used whereby all windows in 
a peak are defined as a cluster. Promoter regions or gene bodies can also be used to cluster windows.

The importance of each window in each cluster can be adjusted by supplying different relative \code{weight} values. 
These weights will be applied during the calculation of the combined p-value. This may be useful for downweighting 
results in low-confidence regions e.g. repeats. Note that these weights are not used to adjust the computed FDR
as weights per cluster are separately defined.
}

\examples{ 
counts<-matrix(rnbinom(400, mu=10, size=20), ncol=4)
y<-DGEList(counts)
design<-model.matrix(~factor(c("a", "a", "b", "b")))
fit<-glmFit(y, design, dispersion=0.05)
out<-glmLRT(fit)
ids<-round(runif(100, 1, 10))
combined<-combineFDR(ids, out$table)
head(combined)

# With window weighting.
w<-round(runif(100, 1, 5))
combined<-combineFDR(ids, out$table, weight=w)
head(combined)

# With cluster weighting, defined as the mean of window weights.
cluster.w<-as.numeric(by(w, list(ids), FUN=mean))
o<-order(combined$PValue)
total.w<-sum(cluster.w)
cs.w<-cumsum(cluster.w[o])
combined$FDR[o]<-rev(cummin(rev(combined$PValue[o]*total.w/cs.w)))
}

\seealso{
\code{\link{mergeWindows}}, \code{\link{glmLRT}}
}

\author{Aaron Lun}

\references{
Simes RJ (1986). An improved Bonferroni procedure for multiple tests of significance. \emph{Biometrika} 73, 751-754.

Benjamini Y and Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society Series B} 57, 289-300. 
}
