\name{checkBimodality}
\alias{checkBimodality}

\title{Check bimodality of regions}

\description{Compute the maximum bimodality score across all base pairs in each region.}

\usage{
checkBimodality(bam.files, regions, width=500, param=readParam(), prior.count=2) 
}

\arguments{
\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
\item{regions}{a \code{GRanges} object specifying the regions over which bimodality is to be calculated}
\item{width}{an integer vector indicating the span with which to compute bimodality}
\item{param}{a \code{readParam} object containing read extraction parameters}
\item{prior.count}{a numeric scalar specifying the prior count to compute bimodality scores}
}

\details{
Consider a base position \code{x}.
This function counts the number of forward- and reverse-strand reads within the interval \code{[x-width+1, x]}.
It then calculates the forward:reverse ratio after adding \code{prior.count} to both counts.
This is repeated for the interval \code{[x, x+width-1]}, and the reverse:forward ratio is then computed.
The smaller of these two ratios is used as the bimodality score.

Sites with high bimodality scores will be enriched for forward- and reverse-strand enrichment on the left and right of the site, respectively.
Given a genomic region, this function will treat each base position as a site.
The largest bimodality score across all positions will be reported for each region.
The idea is to assist with the identification of transcription factor binding sites, which exhibit strong strand bimodality.

If multiple \code{bam.files} are specified, they are effectively pooled so that counting uses all reads in all files.
The value of \code{width} should be set to the average fragment length, as this is allows optimal counting of reads.
Specifically, all reads in the subpeaks on either strand will be counted when \code{x} is set to the center of the peak.
A separate value can be specified for each library, to account for differences in fragmentation.

% Note that the entirety of each read is used for counting, to guarantee symmetry for both forward- and reverse-strand reads.
% This is to avoid high bimodality scores when forward-/reverse-strand reads are stacked directly on top of each other.
% Such read stacks are likely to be mapping artifacts rather than genuine binding sites.
% More generally, symmetry avoids spuriously high scores when you have very high enrichment; if you shift by the read length, there's an opening to get one strand and not the other.

Obviously, this function will be less useful for broad targets like histone marks.
It will also return an error when \code{pe="both"} for \code{param}, as fragments from paired-end data are not stranded.
}

\value{
A numeric vector containing the maximum bimodality score across all bases in each region.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
incoming <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(100, 580, 500, 1500)))

checkBimodality(bamFiles, incoming)
checkBimodality(bamFiles, incoming, width=c(100, 200))
checkBimodality(bamFiles, incoming, param=readParam(minq=20, dedup=TRUE))
checkBimodality(bamFiles, incoming, prior.count=5)
}

\keyword{diagnostics}

