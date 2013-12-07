\name{windowCounts}
\alias{windowCounts}

\title{Count reads overlapping each window}

\description{Count the number of extended reads overlapping a sliding window at spaced positions across the genome.}

\usage{
windowCounts(bam.files, spacing=50, left=0, right=0, ext=100, 
	pet=FALSE, max.frag=500, filter=NULL, bin=NULL, dedup=FALSE, 
	minq=0, restrict=NULL)
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{spacing}{an integer scalar specifying the distance between consecutive windows}
	\item{left}{an integer scalar specifying the width of the window to the left of centre}
	\item{right}{an integer scalar specifying the width of the window to the right of centre}
	\item{ext}{an integer scalar describing the average length of the sequenced fragment}
	\item{pet}{a boolean specifying whether the data involves paired end tags}
	\item{max.frag}{an integer scalar, specifying the maximum allowable fragment size for PET data}
	\item{filter}{an integer scalar for the minimum count sum across libraries for each window}
	\item{bin}{an integer scalar specifying the bin width, if binning is desired}
	\item{dedup}{a boolean specifying whether marked duplicates should be removed}
	\item{minq}{an integer scalar defining the minimum threshold for the mapping quality score}
    \item{restrict}{a character vector specifying which chromosomes should be used for read counting}
}

\value{
A list is returned containing:
	\item{counts}{an integer matrix containing counts for each library at each window remaining after filtering}
    \item{totals}{an integer vector containing the total number of reads in each library}
	\item{region}{a \code{GRanges} object specifying the intervals of each window remaining after filtering}
}

\details{
Each window is defined as the interval \code{[x-left, x+right]} for a centre position \code{x=1, 1+n, 1+2*n, ...} where \code{n=spacing}. 
Each read in \code{bam.files} is extended by \code{ext} in the direction of the read to obtain a rough estimate of the fragment boundaries. 
The number of fragments overlapping the window for each library is then counted for each window position.

The value of \code{ext} can be estimated using \code{correlateReads} or from post-fragmentation gel electrophoresis. Of course, this is
irrelevant when \code{pet=TRUE} where the region bounded by the 5' ends of the reads in each valid pair exactly defines the boundaries 
of the fragment. Fragment sizes above \code{max.frag} are removed. An appropriate choice for this value can be found using \code{getPETSizes}. 

The \code{left} and \code{right} parameters define the width of the binding site i.e. the size of the contact region between the DNA and protein. In 
practical terms, it determines the spatial resolution of the analysis. Larger values for either mean that reads are counted over a larger region 
which results in larger counts. This results in greater detection power at the cost of resolution. 

New windows are found by sliding the current window in the 3' direction by \code{spacing}. Increasing the value of \code{spacing} will reduce 
the number of counts which are extracted. This results in some loss of resolution but it may be necessary when machine memory is limited. 

If \code{filter=NULL}, the count sum filter threshold is defined as the number of libraries multiplied by 5. This removes windows with very low 
counts which cannot be accurately modelled. Increasing \code{filter} can also reduce memory usage by not returning empty or near-empty windows.

If \code{bin!=NULL}, the settings are internally adjusted so that all reads are counted into non-overlapping adjacent bins of width \code{bin}. 
Specifically, \code{right} is set to \code{bin-1}, \code{left} is set to 0, \code{spacing} is set to \code{bin} and \code{filter} is set to 1. 
Only the 5' end of each read or fragment (for paired-end data) is used in counting.

Duplicates must be marked using a tool like Picard's \code{MarkDuplicates} if it is to be removed with \code{dedup=TRUE}. Reads can also be filtered 
by mapping quality score if \code{minq!=0}. However, if \code{pet=TRUE}, only the quality score of the read on the forward strand will be used.  
More generally, users should ensure that PET data is properly synchronised prior to counting. 
}

\author{Aaron Lun}

\seealso{
\code{\link{correlateReads}}, \code{\link{getPETSizes}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles<-system.file("exdata", "chip", c("rep1_sort.bam", "rep2_sort.bam"), package="csaw")
windowCounts(bamFiles, filter=1)
windowCounts(bamFiles, right=100, filter=1)
windowCounts(bamFiles, spacing=100, filter=1)

# Loading PET data.
bamFile<-system.file("exdata", "pet", "pet_sort.bam", package="csaw")
windowCounts(bamFile, pet=TRUE, filter=1)
windowCounts(bamFile, max.frag=100, pet=TRUE, filter=1)
windowCounts(bamFile, max.frag=100, pet=TRUE, filter=1, restrict="chrA")
}
