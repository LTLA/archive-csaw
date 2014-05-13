\name{windowCounts}
\alias{windowCounts}
\alias{countWindows}

\title{Count reads overlapping each window}

\description{Count the number of extended reads overlapping a sliding window at spaced positions across the genome.}

\usage{
windowCounts(bam.files, spacing=50, width=1, ext=100, shift=0,
	pet=c("none", "both", "first", "second"), max.frag=500, rescue.pairs=FALSE,
	filter=NULL, bin=FALSE, dedup=FALSE, minq=0, restrict=NULL, discard=NULL)
countWindows(param, ...)
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{spacing}{an integer scalar specifying the distance between consecutive windows}
	\item{width}{an integer scalar specifying the width of the window}
	\item{ext}{an integer scalar describing the average length of the sequenced fragment}
	\item{shift}{an integer scalar specifying how much the start of each window should be shifted to the left}
	\item{pet}{a character string specifying whether the data involves paired end tags}
	\item{max.frag}{an integer scalar, specifying the maximum allowable fragment size for PET data}
	\item{rescue.pairs}{a logical scalar indicating whether invalid read pairs should be counted}
	\item{filter}{an integer scalar for the minimum count sum across libraries for each window}
	\item{bin}{an integer scalar indicating whether binning should be performed}
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


\section{The sliding window strategy}{
Each window is defined as the interval \code{[x-shift, x+width-1-shift]} for
\code{x=1, 1+spacing, 1+2*spacing, ...}, etc.  The value of \code{width} can be
interpreted as the width of the contact area between the DNA and protein. In
practical terms, it determines the spatial resolution of the analysis. Larger
windows count reads over a larger region which results in larger counts. This 
results in greater detection power at the cost of resolution.

New windows are found by sliding the current window in the 3' direction by
\code{spacing}. Increasing the value of \code{spacing} will reduce the number
of counts which are extracted. This results in some loss of resolution but it
may be necessary when machine memory is limited. The location of each window
can be also tweaked by using \code{shift} to move all windows left (positive values)
or right (negative values).

Each read in \code{bam.files} is extended by \code{ext} in the direction of the
read to obtain a rough estimate of the fragment boundaries.  The number of
fragments overlapping the window for each library is then counted for each
window position. For single-end data, the value of \code{ext} can be estimated
using \code{correlateReads} or from post-fragmentation gel electrophoresis.
For paired-end data, the fragment boundaries can be computed exactly with
\code{pet="both"} for proper pairs.

If \code{bin} is set, settings are internally adjusted so that all reads are
counted into non-overlapping adjacent bins of size \code{width}.  Specifically,
\code{spacing} is set to \code{bin} and \code{filter} is set to 1.  Only the 5'
end of each read or fragment (for paired-end data) is used in counting.
}

\section{Filtering and pruning}{
Windows will be removed if the count sum across all libraries is below
\code{filter}.  This reduces the size of the output by not returning empty or
near-empty windows, which are not usually interesting anyway. If
\code{filter=NULL}, the count sum filter threshold is automatically defined as
the number of libraries multiplied by 5. 

If \code{restrict} is supplied, reads will only be extracted for the specified
chromosomes. This is useful to restrict the analysis to interesting chromosomes
(e.g., no contigs, mitochondria). Conversely, if \code{discard} is set, a read
will be removed if the corresponding alignment is wholly contained within the 
supplied ranges. This is useful for removing reads in repeat regions. 

Marked duplicate reads can be removed with \code{dedup=TRUE}, though the
duplicate field must be properly set beforehand. Reads can also be filtered by
mapping quality score if \code{minq!=0}. This is generally recommended to
remove low-confidence alignments.  The exact threshold for \code{minq} will
depend on the range of scores provided by the aligner. 
}

\section{Additional advice for paired-end data}{
For \code{pet="both"}, users should ensure that the BAM files are properly
synchronised prior to counting. Proper pairs are defined as those reads that
are close together and in an inward-facing orientation. Fragment sizes above
\code{max.frag} are removed; use \code{\link{getPETSizes}} to pick an
appropriate value. Both reads of a pair must also have a mapping quality score
above \code{minq}. 

If \code{rescue=TRUE}, the function will attempt to use proper pairs. If this
is not possible for a read pair, the read with the higher MAPQ score will be
taken, extended to \code{ext} and used as the fragment. This may be useful in
cases where only one read can be mapped from a pair (e.g., if \code{minq} or
\code{discard} is too stringent). Note that both reads will be used for
interchromosomal read pairs.

Finally, paired-end data can also be treated as single-end data by only using
one read from each pair with \code{pet="first"} or \code{"second"}.  This is
useful for poor-quality data where the paired-end procedure has obviously
failed, e.g., many interchromosomal read pairs or pairs with large fragment
lengths.
}

\section{Repeated calls with the same parameters}{
The \code{countWindows} function is a wrapper which allows users to pass in a
list of non-default arguments. This ensures that multiple calls to
\code{windowCounts} can be synchronised easily by passing in the same list in
\code{param}. In particular, users should use the same values for \code{pet},
\code{dedup}, \code{minq}, \code{bam.files} and (for paired-end data)
\code{max.frag}. This ensures that the same reads are being counted in each
call. Note that any additional arguments to \code{countWindows} will override
those specified in \code{param}.
}

\author{Aaron Lun}

\seealso{
\code{\link{correlateReads}}, \code{\link{getPETSizes}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles <- system.file("exdata", "chip", c("rep1_sort.bam", "rep2_sort.bam"), package="csaw")
windowCounts(bamFiles, filter=1)
windowCounts(bamFiles, width=100, filter=1)
windowCounts(bamFiles, spacing=100, filter=1)

# Performing additional non-default calls.
param <- list(bam.files=bamFiles, minq=10, dedup=FALSE)
countWindows(param)
countWindows(param, ext=110, filter=1)
countWindows(param, width=1000, filter=2)
countWindows(param, restrict="chrB", filter=3)
countWindows(param, filter=3, discard=GRanges("chrA", IRanges(1, 100)))

# Loading PET data.
bamFile <- system.file("exdata", "pet", "pet_sort.bam", package="csaw")
windowCounts(bamFile, pet="both", filter=1)
windowCounts(bamFile, pet="first", filter=1)
windowCounts(bamFile, max.frag=100, pet="both", filter=1)
windowCounts(bamFile, max.frag=100, pet="both", filter=1, restrict="chrA")

# Running rescues of PET data (use -1 to coerce single-endedness).
windowCounts(bamFile, max.frag=50, pet="both", filter=1, rescue.pairs=TRUE)
windowCounts(bamFile, max.frag=-1, pet="both", filter=1, rescue.pairs=TRUE)
windowCounts(bamFile, max.frag=100, pet="both", filter=1, restrict="chrA", rescue.pairs=TRUE)
}
