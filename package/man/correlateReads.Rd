\name{correlateReads}
\alias{correlateReads}

\title{Compute correlation coefficients between reads}

\description{Computes the auto- or cross-correlation coefficients between read positions across a set of delay intervals.}

\usage{
correlateReads(bam.files, max.dist=1000, cross=TRUE, param=readParam())
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{max.dist}{integer scalar specifying the maximum delay distance over which correlation coefficients will be calculated}
	\item{cross}{a logical scalar specifying whether cross-correlations should be computed}
	\item{param}{a \code{readParam} object containing read extraction parameters, or a list of such objects (one for each BAM file)}
}

\value{
A numeric vector of length \code{max.dist+1} containing the correlation coefficients for each delay interval from 0 to \code{max.dist}.
}

\details{
If \code{cross=TRUE}, reads are separated into those mapping on the forward and reverse strands. 
Positions on the forward strand are shifted forward by a delay interval. 
The chromosome-wide correlation coefficient between the shifted forward positions and the original reverse positions are computed. 
This is repeated for all delay intervals less than \code{maxDist}. 
A weighted mean for the cross-correlation is taken across all chromosomes, with weighting based on the number of reads. 

Cross-correlation plots can be used to check the quality of immunoprecipitation for ChIP-Seq experiments involving transcription factors or punctate histone marks. 
Strong immunoprecipitation should result in a peak at a delay corresponding to the fragment length. 
A spike may also be observed at the delay corresponding to the read length. 
This is probably an artefact of the mapping process where unique mapping occurs to the same sequence on each strand. 

The construction of cross-correlation plots is usually uninformative for full paired-end data.
This is because the presence of valid pairs will inevitably result in a strong peak at the fragment length.
Nonetheless, immunoprecipitation efficiency can be diagnosed by treating paired-end data as single end data. 
This is done by using only the first or second read based on the value of \code{pe} used in \code{\link{readParam}}.
Setting \code{pe="both"} will result in failure.

If multiple BAM files are specified in \code{bam.files}, the reads from all libraries are pooled prior to calculation of the correlation coefficients.
This is convenient for determining the average correlation profile across an entire dataset. 
Separate calculations for each file will require multiple calls to \code{\link{correlateReads}}.

If \code{cross=FALSE}, auto-correlation coefficients are computed without use of strand information. 
This is designed to guide estimation of the average width of enrichment for diffuse histone marks. 
For example, the width can be defined as the delay distance at which the autocorrelations become negligble.
However, this tends to be ineffective in practice as diffuse marks tend to have very weak correlations to begin with.

By default, marked duplicate reads are removed from each BAM file prior to calculation of coefficients.
This is strongly recommended, even if the rest of the analysis will be performed with duplicates retained.
Otherwise, the read length spike will dominate the plot.
The fragment length peak will no longer be easily visible.
}

\examples{
n <- 20
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
par(mfrow=c(2,2))

x <- correlateReads(bamFile, max.dist=n)
plot(0:n, x, xlab="delay (bp)", ylab="ccf")

x <- correlateReads(bamFile, max.dist=n, param=readParam(dedup=TRUE))
plot(0:n, x, xlab="delay (bp)", ylab="ccf")

x <- correlateReads(bamFile, max.dist=n, cross=FALSE)
plot(0:n, x, xlab="delay (bp)", ylab="acf")
}

\seealso{
	\code{\link{ccf}}
}

\references{
Kharchenko PV, Tolstorukov MY  and Park, PJ (2008). Design and analysis of
ChIP-seq experiments for DNA-binding proteins. \emph{Nat. Biotechnol.} 26,
1351-1359.
}

\author{Aaron Lun}

\keyword{diagnostics}
