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

Note that this does not account for any read extension that might have been performed during read counting.
In such cases, users are advised to expand \code{cur.region} by the extension length on each side.
Counted reads can then be extracted by identifying their extended counterparts that overlap with the original \code{cur.region}.
}

\value{
A \code{GRanges} object is returned.
If \code{pe="both"} in \code{param}, intervals are unstranded and correspond to fragments.
Otherwise, strand-specific intervals that represent reads are returned.
}

\author{Aaron Lun}

\examples{
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
extractReads(GRanges("chrA", IRanges(100, 500)), bamFile)
extractReads(GRanges("chrA", IRanges(100, 500)), bamFile, param=readParam(dedup=TRUE))

bamFile <- system.file("exdata", "pet.bam", package="csaw")
extractReads(GRanges("chrB", IRanges(100, 500)), bamFile)
extractReads(GRanges("chrB", IRanges(100, 500)), bamFile, param=readParam(pe="both"))
extractReads(GRanges("chrB", IRanges(100, 500)), bamFile, param=readParam(pe="first"))

# Dealing with the extension length.
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
ext <- 100
my.reg <- GRanges("chrA", IRanges(10, 200))
my.reg2 <- resize(my.reg, fix="center", width=width(my.reg)+2*ext)
collected <- extractReads(my.reg2, bamFile)

expanded <- resize(collected, width=ext) 
strand(expanded) <- "*"
relevant <- overlapsAny(expanded, my.reg)
collected[relevant,]
}

\seealso{
\code{\link{readParam}}
}

\keyword{visualization}
