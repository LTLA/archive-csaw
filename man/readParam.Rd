\name{readParam}
\docType{class}

\alias{readParam}
\alias{readParam-class}
\alias{$,readParam-method}
\alias{show,readParam-method}
\alias{reform}
\alias{reform,readParam-method}

\title{readParam class and methods}
\description{Class to specify read loading parameters}

\details{
Each readParam object stores a number of parameters, each pertaining to the extraction of reads from a BAM file.
Slots are defined as:
\describe{
	\item{\code{pe}:}{a character string indicating whether paired-end data is present; set to \code{"none"}, \code{"both"}, \code{"first"} or \code{"second"}}
	\item{\code{max.frag}:}{an integer scalar, specifying the maximum fragment length corresponding to a read pair}
	\item{\code{dedup}:}{a logical scalar indicating whether marked duplicate reads should be ignored}
	\item{\code{minq}:}{an integer scalar, specifying the minimum mapping quality score for an aligned read}
	\item{\code{forward}:}{a logical scalar indicating whether only forward reads should be extracted}
	\item{\code{restrict}:}{a character vector containing the names of allowable chromosomes from which reads will be extracted}
	\item{\code{discard}:}{a \code{GRanges} object containing intervals in which any alignments will be discarded}
    \item{\code{BPPARAM}:}{a \code{BiocParallelParam} object specifying if and how parallelization is performed}
}
}

\section{Removing low-quality or irrelevant reads}{
Marked duplicate reads will be removed with \code{dedup=TRUE}.
This may be necessary when many rounds of PCR have been performed during library preparation.
However, it is not recommended for routine counting as it will interfere with the downstream statistical methods.
Note that the duplicate field must be set beforehand in the BAM file for this argument to have any effect.

Reads can also be filtered by their mapping quality scores if \code{minq} is specified at a non-\code{NA} value.  
This is generally recommended to remove low-confidence alignments. 
The exact threshold for \code{minq} will depend on the range of scores provided by the aligner. 
If \code{minq=NA}, no filtering on the score will be performed.

If \code{restrict} is supplied, reads will only be extracted for the specified chromosomes. 
This is useful to restrict the analysis to interesting chromosomes, e.g., no contigs/scaffolds or mitochondria. 
Conversely, if \code{discard} is set, a read will be removed if the corresponding alignment is wholly contained within the supplied ranges. 
This is useful for removing reads in repeat regions. 
}

\section{Parameter settings for paired-end data}{
For \code{pe="both"}, reads are extracted with the previously described filters, i.e., \code{discard}, \code{minq}, \code{dedup}.
Extracted reads are then grouped into proper pairs.
Proper pairs are those where the two reads are close together (on the same chromosome, obviously) and in an inward-facing orientation. 
The fragment interval is defined as that bounded by the 5' ends of the two reads in a proper pair.

The fragment size is defined as the length of the interval bounded by the 5' ends of two inward-facing reads.
Those pairs with fragment sizes above \code{max.frag} are removed, as they are more likely to be the result of mapping errors than genuinely large fragment sizes.
Users should run \code{\link{getPESizes}} to pick an appropriate value for their data sets, though thresholds of around 500-1000 bp are usually fine.

Paired-end data can also be treated as single-end data by only using one read from each pair with \code{pe="first"} or \code{"second"}.  
This is useful for poor-quality data where the paired-end procedure has obviously failed, e.g., with many interchromosomal read pairs or pairs with large fragment lengths.
Treating the data as single-end may allow the analysis to be salvaged.

In all cases, users should ensure that each BAM file containing paired-end data is properly synchronized prior to count loading.
This can be done using standard tools like FixMateInformation from the Picard suite (http://broadinstitute.github.io/picard).
}

\section{Parameter settings for single-end data}{
If \code{pe="none"}, reads are assumed to be single-end.
Read extraction from BAM files is performed with the same quality filters described above.
If \code{forward} is \code{NA}, reads are extracted from all strands.
Otherwise, reads are only extracted from the forward or reverse strands for \code{TRUE} or \code{FALSE}, respectively.
This may be useful for applications requiring strand-specific counting.
A special case is \code{forward=NULL} - see \code{\link{strandedCounts}} for more details.

Any soft clipping in alignments are ignored during extraction (this is also true for paired-end data).
Soft clips are presumed to be sequencing artifacts, e.g., when the adaptor or barcodes are not properly removed from the read sequence.
They are not relevant to computing genomic coverage.
Thus, in this package, any references to the length or 5'/3' ends of the read will actually refer to that of the \emph{alignment}.
This is often more appropriate, e.g., the 5' end of the alignment represents the end of the fragment after clipping of the artifacts.
}

\section{Parallelization options}{
Parallelized read loading can be turned on in several functions by setting \code{BPPARAM} appropriately.
This will usually extract reads from multiple files simultaneously to speed up processing.
Users are referred to \code{?\link{BiocParallelParam}} for more details on how to set \code{BPPARAM}.
By default, a \code{SerialParam} object is used, i.e., no parallelization is performed.
This avoids slow performance on systems that do not support concurrent disk reads.
}

\section{Constructor}{
\describe{
	\item{}{\code{readParam(pe="none", max.frag=500, dedup=FALSE, minq=NA, forward=NA, restrict=NULL, discard=GRanges(), BPPARAM=SerialParam())}:
		Creates a readParam object. 
		Each argument is placed in the corresponding slot, with coercion into the appropriate type.
	}
}
}

\section{Subsetting}{
In the code snippes below, \code{x} is a readParam object.
\describe{
	\item{}{\code{x$name}: 
	Returns the value in slot \code{name}.
	}
}
}

\section{Other methods}{
In the code snippes below, \code{x} is a readParam object.
\describe{
    \item{}{\code{show(x)}:
	Describes the parameter settings in plain English.
  	}
	\item{}{\code{reform(x, ...)}:
		Creates a new readParam object, based on the existing \code{x}.
		Any named arguments in \code{...} are used to modify the values of the slots in the new object, with type coercion as necessary.
	}
}	
}

\author{
Aaron Lun
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}},
\code{\link{correlateReads}},
\code{\link{getPESizes}},
\code{\link{BiocParallelParam}}
}

\examples{
blah <- readParam()
blah <- readParam(discard=GRanges("chrA", IRanges(1, 10)))
blah <- readParam(restrict='chr2')
blah$pe
blah$dedup

# Use 'reform' if only some arguments need to be changed.
blah
reform(blah, dedup=TRUE)
reform(blah, pe="both", max.frag=212.0)
}

\keyword{counting}
