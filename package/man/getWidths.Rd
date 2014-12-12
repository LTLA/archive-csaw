\name{getWidths}
\alias{getWidths}

\title{Get region widths}
\description{Get the widths of the read counting interval for each region.}

\usage{
getWidths(data)
}

\arguments{
\item{data}{a \code{SummarizedExperiment} object, produced by \code{\link{windowCounts}} or \code{\link{regionCounts}}}
}

\details{
Widths of all regions are increased by the average fragment length during the calculations.
This is because each count represents the number of (imputed) fragments overlapping each region.
Thus, a 1 bp window has an effective width that includes the average length of each fragment. 

% It's a direct addition of width and fragment length. If you have a 100 bp interval,
% and you slide a 50 bp fragment along the genome until it hits the interval, there
% are 149 positions where the fragment overlaps the interval. So, the effective width
% is actually 149 bp.

The average fragment length is defined as \code{ext} in the parameter list from \code{paramList(data)}.
If different \code{ext} is used between libraries, the mean value across the dataset is used.
For paired-end data, the average fragment length needs to be specified in \code{ext}.
One recommended strategy is to use the median of the values obtained with \code{\link{getPETSizes}}.
}

\value{
An integer vector containing the effective width, in base pairs, of each region.
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}}
}

\author{
Aaron Lun
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, filter=1)
getWidths(data)

bamFile <- system.file("exdata", "pet.bam", package="csaw")
data <- windowCounts(bamFile, param=readParam(pet="both"), filter=1)
getWidths(data)

# Avoid warning by defining 'ext'.
exptData(data)$param <- reformList(exptData(data)$param, rescue.ext=200)
getWidths(data)
}

\keyword{filtering}
