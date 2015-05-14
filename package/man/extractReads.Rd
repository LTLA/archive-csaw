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

Note that this does not account for any read extension that might have been performed during read counting.
In such cases, users are advised to expand \code{region} by the extension length on each side.
Counted reads can then be extracted by identifying their extended counterparts that overlap with the original \code{region}.

Any strandedness of \code{region} is ignored.
If strand-specific extraction is desired, this can be done by setting \code{param$forward} via \code{\link{reform}}.
Alternatively, the returned \code{GRanges} can be filtered to retain only the desired strand.

If \code{ext} is not \code{NA}, directional read extension will be performed for single-end data.
All extended reads overlapping \code{region} will then be extracted, and the intervals of the extended will be reported..
See \code{\link{windowCounts}} for more details.

By default, unstranded fragments are returned for paired-end data.
If \code{as.reads=TRUE}, the constituent reads are returned, with pairing indicated by \code{pair} in the metadata.
Rescued reads are marked as those with \code{pair} of zero; see \code{\link{readParam}} for details on rescuing.
}

\value{
A \code{GRanges} object is returned.
If \code{pe="both"} in \code{param}, intervals are unstranded and correspond to fragments.
Otherwise, strand-specific intervals that represent reads are returned.
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
    param=readParam(pe="both", rescue.ext=100))
extractReads(bamFile, GRanges("chrB", IRanges(100, 500)), 
    param=readParam(pe="first"))

# Extracting as reads.
extractReads(bamFile, GRanges("chrB", IRanges(100, 500)), 
    param=readParam(pe="both"), as.reads=TRUE)
extractReads(bamFile, GRanges("chrB", IRanges(100, 500)), 
    param=readParam(pe="both", rescue.ext=100), as.reads=TRUE)

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
