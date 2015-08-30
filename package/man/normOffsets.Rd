\name{normOffsets}
\alias{normOffsets}
\alias{normOffsets,matrix-method}

\title{Normalize counts between libraries}

\description{Calculate normalization factors or offsets using count data from multiple libraries.}

\usage{
\S4method{normOffsets}{matrix}(object, lib.sizes=NULL, type=c("scaling", "loess"), 
    weighted=FALSE, ...)
}

\arguments{
  \item{object}{a matrix of integer counts with one column per library}
  \item{lib.sizes}{a numeric vector specifying the total number of reads per library}
  \item{type}{a character string indicating what type of normalization is to be performed}
  \item{weighted}{a logical scalar indicating whether precision weights should be used for TMM normalization}
  \item{...}{other arguments to be passed to \code{\link{calcNormFactors}} for \code{type="scaling"}, or \code{\link{loessFit}} for \code{type="loess"}}
}

\details{
If \code{type="scaling"}, this function provides a convenience wrapper for the \code{\link{calcNormFactors}} function in the \code{edgeR} package.
Specifically, it uses the trimmed mean of M-values (TMM) method to perform normalization. 
Precision weighting is turned off by default so as to avoid upweighting high-abundance regions. 
These are more likely to be bound and thus more likely to be differentially bound. 
Assigning excessive weight to such regions will defeat the purpose of trimming when normalizing the coverage of background regions.

% We can't min/maximize some quantity with regards to the estimated normalization factors, and then choose a reference
% from that, as that would bias the resulting estimates (or at least, interact with existing bias unpredictably).

% Large changes in the normalization factor estimates with changes in the prior suggest that the counts are too low i.e. not
% enough new information in the dataset. This can be overcome by (obviously) increasing the counts. For example, binning
% can be performed with a larger bin size in \code{windowCounts} to obtain proportionally larger counts.

If \code{type="loess"}, this function performs non-linear normalization similar to the fast loess algorithm in \code{\link{normalizeCyclicLoess}}. 
For each sample, a lowess curve is fitted to the log-counts against the log-average count. 
The fitted value for each bin pair is used as the generalized linear model offset for that sample. 
The use of the average count provides more stability than the average log-count when low counts are present for differentially bound regions.

If \code{lib.sizes} is not specified, a warning is issued and the column sums of \code{counts} are used instead.
Note that the same \code{lib.sizes} should be used throughout the analysis if \code{\link{normOffsets}} is called multiple times on the same libraries, e.g., with different bin or window sizes or after different filtering steps.
This ensures that the normalization factors or offsets are comparable between calls.
}

\value{
For \code{type="scaling"}, a numeric vector containing the relative normalization factors for each library.

For \code{type="loess"}, a numeric matrix of the same dimensions as \code{counts}, containing the log-based offsets for use in GLM fitting.
}

\author{Aaron Lun}

\references{
Robinson MD, Oshlack A (2010). A scaling normalization method for differential
expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.

Ballman KV, Grill DE, Oberg AL, Therneau TM (2004). Faster cyclic loess:
normalizing RNA arrays via linear models. \emph{Bioinformatics} 20, 2778-86.
}

\examples{
# A trivial example
counts <- matrix(rnbinom(400, mu=10, size=20), ncol=4)
normOffsets(counts)
normOffsets(counts, lib.sizes=rep(400, 4))

# Adding undersampling
n <- 1000L
mu1 <- rep(10, n)
mu2 <- mu1
mu2[1:100] <- 100
mu2 <- mu2/sum(mu2)*sum(mu1)
counts <- cbind(rnbinom(n, mu=mu1, size=20), rnbinom(n, mu=mu2, size=20))
actual.lib.size <- rep(sum(mu1), 2)
normOffsets(counts, lib.sizes=actual.lib.size)
normOffsets(counts, logratioTrim=0.4, lib.sizes=actual.lib.size)
normOffsets(counts, sumTrim=0.3, lib.size=actual.lib.size)

# With and without weighting, for high-abundance spike-ins.
n <- 100000
blah <- matrix(rnbinom(2*n, mu=10, size=20), ncol=2)
tospike <- 10000
blah[1:tospike,1] <- rnbinom(tospike, mu=1000, size=20)
blah[1:tospike,2] <- rnbinom(tospike, mu=2000, size=20)
full.lib.size <- colSums(blah)

normOffsets(blah, weighted=TRUE, lib.sizes=full.lib.size)
normOffsets(blah, lib.sizes=full.lib.size)
true.value <- colSums(blah[(tospike+1):n,])/colSums(blah)
true.value <- true.value/exp(mean(log(true.value)))
true.value

# Using loess-based normalization, instead.
offsets <- normOffsets(counts, type="loess", lib.size=full.lib.size)
head(offsets)
offsets <- normOffsets(counts, type="loess", span=0.4, lib.size=full.lib.size)
offsets <- normOffsets(counts, type="loess", iterations=1, lib.size=full.lib.size)
}

\seealso{
\code{\link{calcNormFactors}},
\code{\link{loessFit}},
\code{\link{normalizeCyclicLoess}}
}

\keyword{normalization}
