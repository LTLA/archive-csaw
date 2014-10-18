\name{windowCounts}
\alias{windowCounts}

\title{Count reads overlapping each window}

\description{Count the number of extended reads overlapping a sliding window at spaced positions across the genome.}

\usage{
windowCounts(bam.files, spacing=50, width=1, ext=100, shift=0,
	filter=NULL, bin=FALSE, param=readParam())
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{spacing}{an integer scalar specifying the distance between consecutive windows}
	\item{width}{an integer scalar specifying the width of the window}
	\item{ext}{an integer scalar describing the average length of the sequenced fragment}
	\item{shift}{an integer scalar specifying how much the start of each window should be shifted to the left}
	\item{filter}{an integer scalar for the minimum count sum across libraries for each window}
	\item{bin}{an integer scalar indicating whether binning should be performed}
	\item{param}{a \code{readParam} object containing read extraction parameters}
}

\value{
A \code{\link{SummarizedExperiment}} object is returned containing one integer matrix.
Each entry of the matrix contains the count for each library (column) at each window (row).
The coordinates of each window are stored as the \code{rowData}.
The total number of reads in each library are stored as \code{totals} in the \code{colData}.
}

\details{
A window is defined as a genomic interval of size equal to \code{width}. 
The value of \code{width} can be interpreted as the width of the contact area between the DNA and protein. 
In practical terms, it determines the spatial resolution of the analysis. Larger windows count reads over a larger region which results in larger counts. 
This results in greater detection power at the cost of resolution.

By default, the first window on a chromosome starts at base position 1.
This can be shifted to the left by specifying an appropriate value for \code{shift}.
New windows are found by sliding the current window to the right by the specified \code{spacing}.
Increasing \code{spacing} will reduce the frequency at which counts are extracted from the genome. 
This results in some loss of resolution but it may be necessary when machine memory is limited. 

Each read in \code{bam.files} is extended by \code{ext} in the direction of the read to obtain a rough estimate of the fragment boundaries. 
The number of fragments overlapping the window for each library is then counted for each window position. 
For single-end data, the value of \code{ext} can be estimated using \code{correlateReads} or from fragment length diagnostics during library preparation. 
For paired-end data, the fragment boundaries can be computed exactly for proper pairs, by specifying \code{pet="both"} in \code{\link{readParam}}.
If \code{rescue.pairs=TRUE}, improperly paired reads are salvaged by directional extension to \code{ext}.

Windows will be removed if the count sum across all libraries is below \code{filter}. 
This reduces the memory footprint of the output by not returning empty or near-empty windows, which are usually uninteresting anyway. 
If \code{filter=NULL}, the count sum filter threshold is automatically defined as the number of libraries multiplied by 5. 
 
If \code{bin} is set, settings are internally adjusted so that all reads are counted into non-overlapping adjacent bins of size \code{width}. 
Specifically, \code{spacing} is set to \code{bin} and \code{filter} is set to 1. 
Only the 5' end of each read or the midpoint of each fragment (for paired-end data) is used in counting.

}

\author{Aaron Lun}

\seealso{
\code{\link{correlateReads}},
\code{\link{readParam}},
\code{\link{SummarizedExperiment}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
windowCounts(bamFiles, filter=1)
windowCounts(bamFiles, width=100, filter=1)
windowCounts(bamFiles, spacing=100, filter=1)

# Loading PET data.
bamFile <- system.file("exdata", "pet.bam", package="csaw")
windowCounts(bamFile, param=readParam(pet="both"), filter=1)
windowCounts(bamFile, param=readParam(pet="first"), filter=1)
windowCounts(bamFile, param=readParam(max.frag=100, pet="both"), filter=1)
windowCounts(bamFile, param=readParam(max.frag=100, pet="both", restrict="chrA"), filter=1)

# Running rescues of PET data (use -1 to coerce single-endedness ).
windowCounts(bamFile, param=readParam(max.frag=50, pet="both", rescue.pairs=TRUE), filter=1)
}

\keyword{counting}
