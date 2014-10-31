\name{apkb}
\alias{apkb}

\title{Abundance-per-kilobase}
\description{Compute the average abundance-per-kilobase for each region.}

\usage{
apkb(data, prior.count=10, pet.len=NULL, widths=NULL, ...) 
}

\arguments{
\item{data}{a \code{SummarizedExperiment} object produced by \code{\link{windowCounts}} or \code{\link{regionCounts}}}
\item{prior.count}{a numeric scalar specifying the prior count to add}
\item{pet.len}{an integer scalar specifying the average fragment length for paired-end data}
\item{widths}{an integer vector containing the widths of the regions in \code{data}}
\item{...}{other arguments, to be passed to \code{\link{aveLogCPM}}}
}

\details{
This function computes the average abundance-per-kilobase for each region in \code{data}.
Adjustment for the width of each region means that the abundances of regions of differing sizes can be compared.
Note that the adjustment assumes a uniform distribution of fragments throughout each region.
This means that the width of each region is directly proportional to the abundance, such that simple scaling will obtain the abundance-per-kilobase.

Widths of all regions are increased by \code{frag.len} during the calculations.
This is because each count represents the number of (imputed) fragments overlapping each region.
Thus, a 1 bp window has an effective width that includes the average length of each fragment. 
For paired-end data, this length needs to be specified in \code{pet.len}.
Of course, if \code{widths} is specified, this overwrites any internal calculation of the width.

% It's a direct addition of width and fragment length. If you have a 100 bp interval,
% and you slide a 50 bp fragment along the genome until it hits the interval, there
% are 149 positions where the fragment overlaps the interval. So, the effective width
% is actually 149 bp. 

Some subtlety is necessary regarding the treatment of the \code{prior.count}.
Specifically, the prior used in \code{\link{aveLogCPM}} must be increased when the regions are larger.
This ensures that the effective prior is the same after the abundance is scaled down.
Otherwise, the use of the same prior would incorrectly result in a smaller abundance-per-kilobase for larger regions.
The \code{prior.count} value can thus be considered as a prior-count-per-kilobase.
}

\value{
A numeric vector of abundance-per-kilobase values, with one entry for each row of \code{data}.
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}},
\code{\link{aveLogCPM}}
}

\author{
Aaron Lun
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, filter=1)

apkb(data)
apkb(data, prior.count=5)
apkb(data, widths=rep(1000, nrow(data)))
apkb(data, dispersion=0)

bamFile <- system.file("exdata", "pet.bam", package="csaw")
data <- windowCounts(bamFile, param=readParam(pet="both"), filter=1)

try(apkb(data))
apkb(data, pet.len=100)
apkb(data, pet.len=200)
}
