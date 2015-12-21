\name{extractReads}
\alias{extractReads}

\title{Extract reads from a BAM file}

\description{Extract reads from a BAM file with the specified parameter settings.}

\usage{
extractReads(bam.file, region, ext=NA, param=readParam(), as.reads=FALSE)
}

\arguments{
	\item{bam.file}{a character string containing the path to a sorted and indexed BAM file}
	\item{region}{a \code{GRanges} object of length 1 describing the region of interest}
	\item{ext}{an integer scalar specifying the fragment length for directional read extension}
	\item{param}{a \code{readParam} object specifying how reads should be extracted}
	\item{as.reads}{a logical scalar indicating whether reads should be returned instead of fragments for paired-end data}
}	

\details{
This function extracts the reads from a BAM file overlapping a given genomic interval. 
The interpretation of the values in \code{param} is the same as that throughout the package.
The aim is to supply the raw data for visualization, in a manner that maintains consistency with the rest of the analysis.

If \code{pe!="both"} in \code{param}, stranded intervals corresponding to the reads will be reported.
If \code{ext} is not \code{NA}, directional read extension will also be performed -- see \code{\link{windowCounts}} for more details.
If \code{pe="both"}, intervals are unstranded and correspond to fragments from proper pairs.

If \code{as.reads=TRUE} and \code{pe="both"}, the reads in each proper pair are returned directly as a GRangesList of length 2.
The two internal GRanges are of the same length and contain the forward and reverse reads for each proper pair in parallel.
In other words, the nth elements of the first and second GRanges represent the nth proper pair.

Any strandedness of \code{region} is ignored.
If strand-specific extraction is desired, this can be done by setting \code{param$forward} via \code{\link{reform}}.
Alternatively, the returned \code{GRanges} can be filtered to retain only the desired strand.
}

\value{
If \code{pe!="both"} or \code{as.reads=FALSE}, a GRanges object is returned containing the read (for single-end data) or fragment intervals (for paired-end data).
If \code{pe="both"} and \code{as.reads=TRUE}, a GRangesList is returned containing the paired reads -- see Details. 
}

\author{Aaron Lun}

\examples{
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
extractReads(bamFile, GRanges("chrA", IRanges(100, 500)))
extractReads(bamFile, GRanges("chrA", IRanges(100, 500)),
    param=readParam(dedup=TRUE))
extractReads(bamFile, GRanges("chrB", IRanges(100, 500)))

bamFile <- system.file("exdata", "pet.bam", package="csaw")
extractReads(bamFile, GRanges("chrB", IRanges(100, 500)), 
    param=readParam(pe="both"))
extractReads(bamFile, GRanges("chrB", IRanges(100, 500)), 
    param=readParam(pe="first"))

# Extracting as reads.
extractReads(bamFile, GRanges("chrB", IRanges(100, 500)), 
    param=readParam(pe="both"), as.reads=TRUE)

# Dealing with the extension length.
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
my.reg <- GRanges("chrA", IRanges(10, 200))
extractReads(bamFile, my.reg)
extractReads(bamFile, my.reg, ext=100)
}

\seealso{
\code{\link{readParam}},
\code{\link{windowCounts}}
}

\keyword{visualization}
