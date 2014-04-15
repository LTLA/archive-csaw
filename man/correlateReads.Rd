\name{correlateReads}
\alias{correlateReads}

\title{Compute correlation coefficients between reads}

\description{Computes the auto- or cross-correlation coefficients between read positions across a set of delay intervals.}

\usage{
correlateReads(bam.files, max.dist=1000, dedup=FALSE, minq=0, cross=TRUE, restrict=NULL, discard=NULL)
}

\arguments{
	\item{bam.files}{character vectors containing paths to sorted and indexed BAM files}
	\item{max.dist}{integer scalar specifying the maximum delay distance over which correlation coefficients will be calculated}
	\item{dedup}{a logical scalar specifying whether or not marked duplicate reads should be removed}
	\item{minq}{an integer scalar defining the minimum threshold for the mapping quality score}
	\item{cross}{a logical scalar specifying whether cross-correlations should be computed}
	\item{restrict}{a character vector specifying chromosomes for which correlations should be computed}
    \item{discard}{a \code{GRanges} object specifying the regions in which reads are to be ignored}
}

\value{
A numeric vector of length \code{max.dist+1} containing the correlation coefficients for each delay interval from 0 to \code{max.dist}.
}

\details{
If \code{cross=TRUE}, reads are separated into those mapping on the forward and reverse strands. Positions on the forward strand are shifted 
forward by a delay interval and the chromosome-wide correlation coefficients between the shifted forward positions and the original reverse 
positions are computed. This is repeated for all delay intervals less than \code{maxDist}. A weighted mean value for the cross-correlation
is taken across all chromosomes.

Cross-correlation plots are designed to identify the quality of immunoprecipitation for ChIP-Seq experiments involving transcription factors or
punctate histone marks. Strong immunoprecipitation should result in a peak at a delay corresponding to the fragment length. A spike may also be observed 
at the delay corresponding to the read length. This is probably an artefact of the mapping process where unique mapping occurs to the same sequence on each strand. 

If multiple BAM files are specified in \code{bam.files}, the reads from all libraries are pooled prior to calculation of the correlation coefficients.
This is convenient for determining the average correlation profile across an entire dataset. Separate calculations for each file will require multiple 
calls to \code{correlateChIP}.

If \code{cross=FALSE}, auto-correlation coefficients are computed without use of strand information. This is designed to guide estimation of the 
average width of enrichment for diffuse histone marks. For example, the width can be defined as the delay distance at which the autocorrelations become negligble.
However, this tends to be ineffective in practice as diffuse marks tend to have very weak correlations to begin with.

Note that any duplicate read must be marked in the bit field of the BAM file using a tool like Picard if it is to be removed with \code{dedup=TRUE}.
}

\examples{
n<-20
bamFile<-system.file("exdata", "chip", "rep1_sort.bam", package="csaw")
par(mfrow=c(2,2))

x<-correlateReads(bamFile, max.dist=n)
plot(0:n, x, xlab="delay (bp)", ylab="ccf")

x<-correlateReads(bamFile, max.dist=n, discard=GRanges("chrB", IRanges(1, 100)))
plot(0:n, x, xlab="delay (bp)", ylab="ccf")

x<-correlateReads(bamFile, max.dist=n, restrict="chrA")
plot(0:n, x, xlab="delay (bp)", ylab="ccf")

x<-correlateReads(bamFile, max.dist=n, cross=FALSE)
plot(0:n, x, xlab="delay (bp)", ylab="acf")
}

\seealso{
	\code{\link{ccf}}
}

\references{
Kharchenko PV, Tolstorukov MY  and Park, PJ (2008). Design and analysis of ChIP-seq experiments for DNA-binding proteins. \emph{Nature Biotechnology} 26, 1351-1359.
}

\author{Aaron Lun}
