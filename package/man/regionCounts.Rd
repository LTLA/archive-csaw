\name{regionCounts}
\alias{regionCounts}

\title{Count reads overlapping each region}

\description{Count the number of extended reads overlapping pre-specified regions}

\usage{
regionCounts(bam.files, regions, ext=100, param=readParam())
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{regions}{a \code{GRanges} object containing the regions over which reads are to be counted}
	\item{ext}{an integer scalar describing the average length of the sequenced fragment}
	\item{param}{a \code{readParam} object containing read extraction parameters}
}

\value{
A \code{\link{SummarizedExperiment}} object is returned containing one integer matrix.
Each entry of the matrix contains the count for each library (column) at each region (row).
The coordinates of each region are stored as the \code{rowData}.
The total number of reads in each library are stored as \code{totals} in the \code{colData}.
}

\details{
This function simply provides a wrapper around \code{\link{countOverlaps}} for read counting into specified regions. 
It is provided so as to allow for counting with awareness of the other parameters, e.g., \code{ext}, \code{pet}.
This allows users to coordinate region-based counts with those from \code{\link{windowCounts}}. 
Checking that the output \code{totals} are the same between the two calls is strongly recommended.
}

\author{Aaron Lun}

\seealso{
\code{\link{countOverlaps}},
\code{\link{windowCounts}},
\code{\link{readParam}},
\code{\link{SummarizedExperiment}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
incoming <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(200, 1000, 700, 1500)))
regionCounts(bamFiles, regions=incoming)
regionCounts(bamFiles, regions=incoming, param=readParam(restrict="chrB"))

# Loading PET data.
bamFile <- system.file("exdata", "pet.bam", package="csaw")
regionCounts(bamFile, regions=incoming, param=readParam(pet="both"))
regionCounts(bamFile, regions=incoming, param=readParam(max.frag=100, 
	pet="first", restrict="chrA"))
regionCounts(bamFile, regions=incoming, param=readParam(max.frag=100, 
	pet="both", restrict="chrA", rescue.pairs=TRUE))
}

\keyword{counting}
