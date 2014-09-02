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
Each readParam object stores a number of parameters to extract reads from a BAM file.
Slots are defined as:
\describe{
	\item{\code{pet}:}{a character string indicating whether paired-end data is present}
	\item{\code{max.frag}:}{an integer scalar, specifying the maximum fragment length corresponding to a read pair}
	\item{\code{rescue.pairs}:}{a logical scalar indicating whether invalid read pairs should be rescued}
	\item{\code{rescue.ext}:}{an integer scalar indicating the extension length for invalid read pairs}
	\item{\code{dedup}:}{a logical scalar indicating whether marked duplicate reads should be ignored}
	\item{\code{minq}:}{an integer scalar, specifying the minimum mapping quality score for an aligned read}
	\item{\code{restrict}:}{a character vector containing the names of allowable chromosomes from which reads will be extracted}
	\item{\code{discard}:}{a \code{GRanges} object containing intervals in which any alignments will be discarded}
}

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

The \code{pet} parameter can take values of \code{"none"}, \code{"both"}, \code{"first"} or \code{"second"}.
If \code{pet="none"}, reads are assumed to be single-end.
Each read will then be processed on its own merits.
The other settings for \code{pet} and the values of \code{max.frag} and \code{rescue.pairs} describe the treatment of paired-end data; see below for more information.
}

\section{Additional advice for paired-end data}{
Reads are first extracted in a single-end manner, using the parameters described above, e.g., \code{discard}, \code{minq}, \code{dedup}.
For \code{pet="both"}, the extracted reads are grouped into proper pairs.
Proper pairs are defined as reads that are close together and in an inward-facing orientation. 
The fragment interval is defined as that bounded by the 5' ends of the two reads in a proper pair.
Fragment sizes above \code{max.frag} are removed; use \code{\link{getPETSizes}} to pick an appropriate value. 

If \code{rescue.pairs=FALSE}, only reads in proper pairs are used to construct a fragment interval.
Otherwise, if \code{rescue.pairs=TRUE}, the function will first attempt to assign reads into proper pairs.
For the remaining reads in improper pairs, the read with the higher MAPQ score will be taken and directionally extended to a length of \code{rescue.ext}.
The extended read will then be used as the fragment interval.
This may be useful in cases where only one read can be mapped from a pair (e.g., if \code{minq} or \code{discard} is too stringent). 
Note that both reads will be extended for interchromosomal read pairs.

Finally, paired-end data can also be treated as single-end data by only using one read from each pair with \code{pet="first"} or \code{"second"}.  
This is useful for poor-quality data where the paired-end procedure has obviously failed, e.g., many interchromosomal read pairs or pairs with large fragment lengths.
Treating the data as single-end may allow the analysis to be salvaged.

In all cases, users should ensure that each BAM file containing paired-end data is properly synchronized prior to count loading.
}

\section{Constructor}{
\describe{
	\item{}{\code{readParam(pet="none", max.frag=500, rescue.pairs=FALSE, rescue.ext=200, dedup=FALSE, minq=NA, restrict=NULL, discard=GRanges())}:
		Creates a readParam object. 
		Each argument is placed in the corresponding slot, with coercion into the appropriate type.
	}
}
}

\section{Subsetting}{
In the code snippets below, \code{x} is a readParam object.
\describe{
	\item{}{\code{x$name}: 
	Returns the value in slot \code{name}.
	}
}
}

\section{Other methods}{
In the code snippets below, \code{x} is a readParam object.
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
\code{\link{extractReads}},
\code{\link{getPETSizes}}
}

\examples{
blah <- readParam()
blah <- readParam(discard=GRanges("chrA", IRanges(1, 10)))
blah <- readParam(restrict='chr2')
blah$pet
blah$dedup

# Use 'reform' if only some arguments need to be changed.
blah
reform(blah, dedup=TRUE)
reform(blah, pet="both", max.frag=212.0)
}

\keyword{counting}
