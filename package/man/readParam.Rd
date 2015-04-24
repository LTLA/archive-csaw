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
	\item{\code{rescue.ext}:}{an integer scalar indicating the extension length for rescued reads from invalid pairs}
	\item{\code{fast.pe}:}{a logical scalar specifying whether fast fragment extraction should be performed for paired-end data}
	\item{\code{dedup}:}{a logical scalar indicating whether marked duplicate reads should be ignored}
	\item{\code{minq}:}{an integer scalar, specifying the minimum mapping quality score for an aligned read}
	\item{\code{forward}:}{a logical scalar indicating whether only forward reads should be extracted}
	\item{\code{restrict}:}{a character vector containing the names of allowable chromosomes from which reads will be extracted}
	\item{\code{discard}:}{a \code{GRanges} object containing intervals in which any alignments will be discarded}
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
Proper pairs are those where the two reads are close together and in an inward-facing orientation. 
The fragment interval is defined as that bounded by the 5' ends of the two reads in a proper pair.
Fragment sizes above \code{max.frag} are removed; use \code{\link{getPESizes}} to pick an appropriate value. 

By default, only reads in proper pairs are used to construct a fragment interval.
If \code{rescue.ext!=NA}, the function will also attempt to recover reads from improper pairs.
For each improper pair, the read with the higher MAPQ score will be directionally extended with \code{rescue.ext}.
Users should set \code{rescue.ext} as the median fragment length from \code{\link{getPESizes}}.
The extended read will then be used as the fragment interval.
Similarly, any reads without a mapped mate will be extended. 
For reasons of convenience, both reads will be rescued and extended for interchromosomal read pairs.
This is acceptable as these reads will never be counted twice into the same region.

Each run through the position-sorted BAM file tends to be time-consuming, as read names must be extracted and matched.
This can be avoided by setting \code{fast.pe=TRUE}.
In this case, only cursory checks are done with regards to proper pairing.
Each fragment is defined from a forward-stranded read with a reverse-stranded mate and a positive insert size no larger than \code{max.frag}.
Values of \code{minq} and \code{discard} are ignored, as the specifics of the mate alignment are not accessed in this fast mode.

Finally, paired-end data can also be treated as single-end data by only using one read from each pair with \code{pe="first"} or \code{"second"}.  
This is useful for poor-quality data where the paired-end procedure has obviously failed, e.g., with many interchromosomal read pairs or pairs with large fragment lengths.
Treating the data as single-end may allow the analysis to be salvaged.

In all cases, users should ensure that each BAM file containing paired-end data is properly synchronized prior to count loading.
}

\section{Parameter settings for single-end data}{
If \code{pe="none"}, reads are assumed to be single-end.
Read extraction from BAM files is performed with the same quality filters described above.
If \code{forward=NA}, reads are extracted from all strands.
Otherwise, reads are only extracted from the forward or reverse strands for \code{TRUE} or \code{FALSE}, respectively.
This may be useful for applications requiring strand-specific counting.
A special case is \code{forward=NULL} - see \code{\link{strandedCounts}} for more details.

Note that directional extension is often used in analyses of single-end data.
However, it must be stressed that the \code{rescue.ext} parameter here does \emph{not} determine the length of the extension.
Instead, the average fragment length must be specified in each of the relevant downstream functions.
This is because directional extension of single-end data pertains to analysis, not read extraction.

% Doing so for paired-end data is just to recover the paired-end fragment, for similarity with proper pairs.
}

\section{Constructor}{
\describe{
	\item{}{\code{readParam(pe="none", max.frag=500, rescue.ext=NA, dedup=FALSE, minq=NA, forward=NA, restrict=NULL, discard=GRanges())}:
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
\code{\link{getPESizes}}
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
reform(blah, rescue.ext=200)
reform(blah, pe="both", max.frag=212.0)
reform(blah, pe="both", fast.pe=TRUE)
}

\keyword{counting}
