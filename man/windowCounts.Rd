\name{windowCounts}
\alias{windowCounts}
\alias{countWindows}

\title{Count reads overlapping each window}

\description{Count the number of extended reads overlapping a sliding window at spaced positions across the genome.}

\usage{
windowCounts(bam.files, spacing=50, left=0, right=0, ext=100, 
	pet=c("none", "both", "first", "second"), max.frag=500, filter=NULL, bin=NULL, dedup=FALSE, 
	minq=0, restrict=NULL, discard=NULL)
countWindows(param, ...)
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{spacing}{an integer scalar specifying the distance between consecutive windows}
	\item{left}{an integer scalar specifying the width of the window to the left of centre}
	\item{right}{an integer scalar specifying the width of the window to the right of centre}
	\item{ext}{an integer scalar describing the average length of the sequenced fragment}
	\item{pet}{a character string specifying whether the data involves paired end tags}
	\item{max.frag}{an integer scalar, specifying the maximum allowable fragment size for PET data}
	\item{filter}{an integer scalar for the minimum count sum across libraries for each window}
	\item{bin}{an integer scalar specifying the bin width, if binning is desired}
	\item{dedup}{a boolean specifying whether marked duplicates should be removed}
	\item{minq}{an integer scalar defining the minimum threshold for the mapping quality score}
    \item{restrict}{a character vector specifying which chromosomes should be used for read counting}
	\item{discard}{a \code{GRanges} object specifying the regions in which reads are to be ignored}
	\item{param}{a named list of non-default arguments with appropriate values}
	\item{...}{other arguments to be passed to \code{windowCounts}}
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

For single-end data, the value of \code{ext} can be estimated using \code{correlateReads} or from post-fragmentation gel electrophoresis. 
For paired-end data, the fragment boundaries can be computed exactly with \code{pet="both"} for valid pairs with reads in an inward-facing 
orientation. Fragment sizes above \code{max.frag} are removed; use \code{getPETSizes} to pick an appropriate value. Poor-quality paired-end
data can be treated as single-end data by only using one read from each pair with \code{pet="first"} or \code{"second"}.

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

If \code{restrict} is set, reads will only be extracted from the chromosomes contained in the supplied vector. This is useful to restrict the
analysis to interesting chromosomes (e.g., no contigs). Conversely, if \code{discard} is set, reads that overlap the supplied ranges are 
discarded. This is useful for removing reads in repeat regions.

The duplicate flag must be set if reads are to be removed with \code{dedup=TRUE}. For \code{pet="both"}, users should ensure that the BAM files 
are properly synchronised prior to counting. Both these procedures can be performed using tools from the Picard suite prior to read counting.	  
	
Reads can also be filtered by mapping quality score if \code{minq!=0}. This is generally recommended to remove low-confidence alignments.
The exact threshold for \code{minq} will depend on the range of scores provided by the aligner. If \code{pet="both"}, both reads of a pair 
must have a mapping quality score above \code{minq}.

The \code{countWindows} function is a wrapper which allows users to pass in a list of non-default arguments. This ensures that multiple calls
to \code{windowCounts} can be synchronised easily by passing in the same list in \code{param}. In particular, users should use the same values
for \code{pet}, \code{dedup}, \code{minq}, \code{bam.files} and (for paired-end data) \code{max.frag}. This ensures that the same reads are 
being counted in each call. Note that any additional arguments to \code{countWindows} will override those specified in \code{param}.
}

\author{Aaron Lun}

\seealso{
\code{\link{correlateReads}}, \code{\link{getPETSizes}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles <- system.file("exdata", "chip", c("rep1_sort.bam", "rep2_sort.bam"), package="csaw")
windowCounts(bamFiles, filter=1)
windowCounts(bamFiles, right=100, filter=1)
windowCounts(bamFiles, spacing=100, filter=1)

# Performing additional non-default calls.
param <- list(bam.files=bamFiles, minq=10, dedup=FALSE)
countWindows(param)
countWindows(param, ext=110, filter=1)
countWindows(param, right=999, filter=2)
countWindows(param, restrict="chrB", filter=3)
countWindows(param, filter=3, discard=GRanges("chrA", IRanges(1, 100)))

# Loading PET data.
bamFile <- system.file("exdata", "pet", "pet_sort.bam", package="csaw")
windowCounts(bamFile, pet="both", filter=1)
windowCounts(bamFile, pet="first", filter=1)
windowCounts(bamFile, max.frag=100, pet="both", filter=1)
windowCounts(bamFile, max.frag=100, pet="both", filter=1, restrict="chrA")
}
