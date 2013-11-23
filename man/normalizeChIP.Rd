\name{normalizeChIP}
\alias{normalizeChIP}

\title{Normalize for undersampling in ChIP-seq data}

\description{Calculate normalization factors to scale the raw library sizes using a modified TMM implementation.}

\usage{
normalizeChIP(counts, lib.sizes=colSums(counts), weighted=FALSE, ...)
}

\arguments{
  \item{counts}{a matrix of integer counts with one column per library}
  \item{lib.sizes}{a numeric vector specifying the total number of reads per library}
  \item{weighted}{a logical scalar indicating whether precision weights should be used for TMM}
  \item{...}{other arguments to be passed to \code{calcNormFactors}}
}

\details{
This function provides a convenience wrapper for the trimmed mean of M-values (TMM) method implemented in the \code{calcNormFactors}
function from the \code{edgeR} package. Briefly, log-expression ratios (M-values) are computed between each library and 
a selected reference library. Trimming is performed to remove M-values which are themselves extreme or correspond to 
extreme abundances. The normalization factor is then defined from the mean of the remaining M-values.

% This implementation differs from the original by providing support for prior count values. The addition of a \code{prior.count}
% prevents undefined M-values when the count for either library is zero. Simply discarding these instances can lead to some 
% unpredictable results, espeically in ChIP-seq data where a non-negligble proportion of counts can be zero. Note that instances 
% involving two zero counts are still discarded.

% Compare the 'low count' example to calcNormFactors, for potential issues with the precision of the estimates.

Note that the precision weighting in the original TMM procedure is turned off by default. The weighting scheme will increase the 
contribution of points where the counts are high. However, differential binding in a ChIP-seq experiment is more likely to occur 
at high abundance features corresponding to enriched regions. If any differentially bound regions survive trimming (e.g. those 
with subtle fold changes), upweighting them will defeat the purpose of trimming to remove potential DE regions in the first place.

% If \code{refColumn=NULL}, the reference library is defined as that with the median count closest to the mean of the median counts
% across all librarites. This is designed to select the library with the most intermediate counts. Such a library will generally be 
% most suited for comparison to all other libraries.

% We can't min/maximize some quantity with regards to the estimated normalization factors, and then choose a reference
% from that, as that would bias the resulting estimates (or at least, interact with existing bias unpredictably).

% Large changes in the normalization factor estimates with changes in the prior suggest that the counts are too low i.e. not
% enough new information in the dataset. This can be overcome by (obviously) increasing the counts. For example, binning
% can be performed with a larger bin size in \code{windowCounts} to obtain proportionally larger counts.
}

\value{
A numeric vector containing the relative normalization factors for each library.
}

\author{Aaron Lun}

\references{
Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
}

\examples{
# A trivial example
counts<-matrix(rnbinom(400, mu=10, size=20), ncol=4)
normalizeChIP(counts)
normalizeChIP(counts, lib.sizes=rep(400, 4))

# Adding undersampling
n<-1000L
mu1<-rep(10, n)
mu2<-mu1
mu2[1:100]<-100
mu2<-mu2/sum(mu2)*sum(mu1)
counts<-cbind(rnbinom(n, mu=mu1, size=20), rnbinom(n, mu=mu2, size=20))
normalizeChIP(counts)
normalizeChIP(counts, logratioTrim=0.4)
normalizeChIP(counts, sumTrim=0.3)

# An example with and without weighting, in the presence of high-abundance spike-ins.
n<-100000
blah<-matrix(rnbinom(2*n, mu=10, size=20), ncol=2)
tospike<-10000
blah[1:tospike,1]<-rnbinom(tospike, mu=1000, size=20)
blah[1:tospike,2]<-rnbinom(tospike, mu=2000, size=20)

normalizeChIP(blah, weighted=TRUE)
normalizeChIP(blah)
true.value<-colSums(blah[(tospike+1):n,])/colSums(blah)
true.value<-true.value/exp(mean(log(true.value)))
true.value
}


\seealso{
\code{\link{calcNormFactors}}
}
