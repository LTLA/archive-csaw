\name{windowCounts}
\alias{windowCounts}

\title{Count reads overlapping each window}

\description{Count the number of extended reads overlapping a sliding window at spaced positions across the genome.}

\usage{
windowCounts(bam.files, spacing=50, width=spacing, ext=100, shift=0,
	filter=10, bin=FALSE, param=readParam())
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{spacing}{an integer scalar specifying the distance between consecutive windows}
	\item{width}{an integer scalar specifying the width of the window}
	\item{ext}{an integer scalar or vector, or a list of two integer scalars/vectors, containing the average length(s) of the sequenced fragments in each library}
	\item{shift}{an integer scalar specifying how much the start of each window should be shifted to the left}
	\item{filter}{an integer scalar for the minimum count sum across libraries for each window}
	\item{bin}{a logical scalar indicating whether binning should be performed}
	\item{param}{a \code{readParam} object containing read extraction parameters, or a list of such objects (one for each BAM file)}
}

\value{
A \code{RangedSummarizedExperiment} object is returned containing one integer matrix. 
Each entry of the matrix contains the count for each library (column) at each window (row).
The coordinates of each window are stored as the \code{rowRanges}.
The total number of reads in each library are stored as \code{totals} in the \code{colData}, along with the read extension length and \code{param} used in each library.
Other window counting parameters (e.g., \code{spacing}, \code{width}) are stored in the \code{metadata}.
}

\section{Defining the sliding windows}{
A window is defined as a genomic interval of size equal to \code{width}. 
The value of \code{width} can be interpreted as the width of the contact area between the DNA and protein. 
In practical terms, it determines the spatial resolution of the analysis. Larger windows count reads over a larger region which results in larger counts. 
This results in greater detection power at the cost of resolution.

By default, the first window on a chromosome starts at base position 1.
This can be shifted to the left by specifying an appropriate value for \code{shift}.
New windows are found by sliding the current window to the right by the specified \code{spacing}.
Increasing \code{spacing} will reduce the frequency at which counts are extracted from the genome. 
This results in some loss of resolution but it may be necessary when machine memory is limited. 

If \code{bin} is set, settings are internally adjusted so that all reads are counted into non-overlapping adjacent bins of size \code{width}. 
Specifically, \code{spacing} is set to \code{width} and \code{filter} is capped at a maximum value of 1 (empty bins can be retained with \code{filter=0}).
Only the 5' end of each read or the midpoint of each fragment (for paired-end data) is used in counting.
}

\section{Read extraction and counting}{ 
Read extraction from the BAM files is governed by the \code{param} argument.
Specifying a single \code{readParam} object will use the same extraction parameters for all files.
Different parameters can also be used for each file by specifying a list of \code{readParam} objects.
However, users should take care with library-specific parameters, lest spurious differences be introduced between libraries.

Fragments are inferred from reads by directional extension (single-end; see below) or by identifying proper pairs (paired-end; see \code{\link{readParam}} for more details).
The number of fragments overlapping the window for each library is then counted for each window position. 
Windows will be removed if the count sum across all libraries is below \code{filter}. 
This reduces the memory footprint of the output by not returning empty or near-empty windows, which are usually uninteresting anyway. 
 
The strandedness of the output \code{rowRanges} is set based on the strand(s) from which the reads are extracted and counted.
This is determined by the value of the \code{forward} slot in the input \code{param} object.
}

\section{Elaborating on directional extension}{
For single-end reads, directional extension is performed whereby each read is extended from its 3' end to the average fragment length, i.e., \code{ext}.
This obtains a rough estimate of the interval of the fragment from which the read was derived.
It is particularly useful for TF data, where extension specifically increases the coverage of peaks that exhibit strand bimodality.
Substantially different fragment lengths between libraries can be accommodated by supplying a vector to \code{ext}, where each entry represents the extension length for the corresponding library.
If any value of \code{ext} is set to \code{NA}, the read length is used as the fragment length in that library.

However, different fragment lengths will result in different peak widths between libraries.
This may result in the detection of irrelevant differences corresponding to these differences in widths.
To avoid this, fragment lengths in all libraries can be scaled to \code{final.ext}.
For a bimodal peak, scaling effectively aligns the subpeaks on a given strand across all libraries to a common location.
This removes the most obvious differences in widths.

The value of \code{final.ext} can be specified by passing a list (or equivalently subsettable object) of integer vectors in the \code{ext} argument.
The library-specific fragment lengths are taken from the first element of the list, while the final fragment length is defined from the second element. 
To this end, the second element should be either a scalar or a vector containing a single duplicated value.
If \code{final.ext} is set to \code{NA}, no rescaling is performed from the library-specific fragment lengths.
}

\section{Comments on ext for paired-end data}{
Directional extension is not performed for paired-end data, so the values in \code{ext} are not used directly.
Hwoever, rescaling can still be performed to standardize fragment lengths across libraries, by resizing each fragment from its midpoint.
This will use the second element of \code{ext} as \code{final.ext}, if \code{ext} is specified as a list of length 2.

On a similar note, some downstream functions will use the extension length in the output \code{colData} as the average fragment length.
Thus, to maintain compatibility, the extension length in the output \code{colData} is set to the average of the inferred fragment lengths for valid read pairs. 
These values will not be used in \code{windowCounts}, but instead, in functions like \code{\link{getWidths}}.
}

\author{Aaron Lun}

\seealso{
\code{\link{correlateReads}},
\code{\link{readParam}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
windowCounts(bamFiles, filter=1)
windowCounts(bamFiles, width=100, filter=1)
windowCounts(bamFiles, ext=c(50, 100), spacing=100, filter=1)
windowCounts(bamFiles, ext=DataFrame(c(50, 100), 80), width=100)

# Loading PE data.
bamFile <- system.file("exdata", "pet.bam", package="csaw")
windowCounts(bamFile, param=readParam(pe="both"), filter=1)
windowCounts(bamFile, param=readParam(pe="first"), filter=1)
windowCounts(bamFile, param=readParam(max.frag=100, pe="both"), filter=1)
windowCounts(bamFile, param=readParam(max.frag=100, pe="both", restrict="chrA"), filter=1)
}

\keyword{counting}
