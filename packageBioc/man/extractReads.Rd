\name{extractReads}
\alias{extractReads}

\title{Extract reads from a BAM file}

\description{Extract reads from a BAM file with the specified parameter settings.}

\usage{
extractReads(cur.region, bam.file, param=readParam())
}

\arguments{
	\item{cur.region}{a \code{GRanges} object of length 1 describing the region of interest}
	\item{bam.file}{a character string containing the path to a sorted and indexed BAM file}
	\item{param}{a \code{readParam} object specifying how reads should be extracted}
}	

\details{
This function extracts the reads from a BAM file overlapping a given genomic interval. 
The interpretation of the values in \code{param} is the same as that throughout the package.
The aim is to supply the raw data for visualization, in a manner that maintains consistency with the rest of the analysis.
}

\value{
A \code{GRanges} object is returned.
If \code{pet="both"} in \code{param}, intervals are unstranded and correspond to fragments.
Otherwise, strand-aware reads are returned.
}

\author{Aaron Lun}

\examples{
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
extractReads(GRanges("chrA", IRanges(100, 500)), bamFile)
extractReads(GRanges("chrA", IRanges(100, 500)), bamFile, param=readParam(dedup=TRUE))

bamFile <- system.file("exdata", "pet.bam", package="csaw")
extractReads(GRanges("chrB", IRanges(100, 500)), bamFile)
extractReads(GRanges("chrB", IRanges(100, 500)), bamFile, param=readParam(pet="both"))
extractReads(GRanges("chrB", IRanges(100, 500)), bamFile, param=readParam(pet="first"))
}

\seealso{
\code{\link{readParam}}
}

\keyword{visualization}
