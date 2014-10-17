\name{profileSummit}
\alias{profileSummit}

\title{Profile binding summits}
\description{Get the coverage profile around the summits of potential binding sites.}

\usage{
profileSummit(bam.files, ext=100, width=5000, res=50, min.depth=1, param=readParam()) 
}

\arguments{
\item{bam.files}{a character vector containing paths to one or more BAM files}
\item{ext}{an integer scalar specifying the average fragment length for single-end data}
\item{width}{an integer scalar specifying the half-window over which the local maxima will be computed}
\item{res}{an integer scalar specifying the window size for coverage summarization}
\item{min.depth}{a numeric scalar defining the minimum number of reads-per-base at a potential summit}
\item{param}{a \code{readParam} object containing read extraction parameters}
}

\details{
This function summarizes the coverage profile around summits of putative binding sites.
The profile can then be used to guide an intelligent choice of window size.
The algorithm itself is a stripped-down version of that used in the MACS peak caller.

The first step involves summarizing the coverage into non-overlapping windows of size \code{res}. 
Single-end reads are extended by \code{ext} in the direction of the read to impute the fragment.
For paired-end reads, the interval between each pair is used as the fragment.
The number of fragments overlapping each window is counted.
If multiple \code{bam.files} are specified, reads are pooled across files for counting into each window.

The next step involves identifying the local maxima. 
For each window, the set of windows within a distance of \code{width} is considered on both sides.
If the center window has a read count above all of its neighbouring windows, it is treated as the local maximum (no ties allowed).
In addition, a local maximum window is only defined as a summit if the count for that window is above \code{min.depth*res}.
This ensures that only high-abundance maxima are considered as potential windowding sites.
The per-base definition of \code{min.depth} ensures that it is consistent if \code{res} is changed.

For each summit, the coverage is recalculated after expanding each the size of the central window.
This is repeated over several expansions, and the coverage at each expanded window size is recorded.
The relative coverage is defined by dividing each coverage by the original (i.e., unexpanded) coverage at the summit.
This is averaged across all summits to obtain the binding profile for the library.
}

\value{
A list of two components:
\item{span}{an integer vector containing the expanded window sizes}
\item{coverage}{a numeric vector of average relative coverages for each expanded window size}
}

\author{
Aaron Lun
}

\examples{
# Reduce min.depth for this example.
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
x <- profileSummit(bamFile, min.depth=0.1, width=1000)
plot(x$span, x$coverage)

x <- profileSummit(bamFile, min.depth=0.1, width=500, param=readParam(dedup=TRUE))
plot(x$span, x$coverage)

x <- profileSummit(bamFile, min.depth=0.1, res=20, width=100)
plot(x$span, x$coverage)
}

\references{
Zhang Y, Liu T, Meyer CA et al. (2008). Model-based Analysis of ChIP-Seq (MACS). \emph{Genome Biol} 9, R137
}

\keyword{diagnostics}
