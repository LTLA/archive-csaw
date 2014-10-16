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
\item{res}{an integer scalar specifying the bin size for coverage summarization}
\item{min.depth}{a numeric scalar defining the minimum number of reads-per-base at a potential summit}
\item{param}{a \code{readParam} object containing read extraction parameters}
}

\details{
This function summarizes the coverage profile around summits of putative binding sites.
The profile can then be used to guide an intelligent choice of window size.
The algorithm itself is a stripped-down version of that used in the MACS peak caller.

The first step involves summarizing the coverage into bins of size \code{res}. 
Single-end reads are shifted by half of \code{ext} in the direction of the read, and the 5' end is used for counting.
For paired-end reads, the midpoint of the fragment associated with each pair is used.
If multiple \code{bam.files} are specified, reads are pooled across files for counting into each bin.

The next step involves identifying the local maxima. 
For each bin, the set of bins within \code{width} is considered on both sides.
If the center bin has coverage above those other bins, it is treated as the local maximum.
Ties are not allowed.

A local maximum bin is only defined as a summit if the count for that bin is above \code{min.depth*res}.
This means that only high-abundance maxima are considered as potential binding sites.
The per-base definition of \code{min.depth} ensures that it is consistent if \code{res} is changed.

For each summit, the relative coverage of a surrounding bin is defined as the ratio of the count for that bin to the count at the summit.
The profile is defined by recording the relative coverage at each distance away from the summit.
This is averaged over both sides of the summit and across all summits.

Users can plot the profile to identify the distance from the summit at which the relative coverage becomes negligble.
The window size can be roughly defined by doubling the chosen distance.
This considers the reads on both sides of each summit during counting.
}

\value{
A list of two components:
\item{span}{an integer vector containing the distances of each surrounding bin from the summit}
\item{coverage}{a numeric vector of average relative coverages for each distance from the summit} 
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
