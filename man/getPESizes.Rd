\name{getPESizes}
\alias{getPESizes}

\title{Compute fragment lengths for paired-end tags}
 
\description{Compute the length of the sequenced fragment for each read pair in paired-end tag (PE) data.}

\usage{
getPESizes(bam.file, param=readParam(pe="both"))
}

\arguments{
	\item{bam.file}{a character string containing the file path to a sorted and indexed BAM file}
    \item{param}{a \code{readParam} object containing read extraction parameters}
}

\value{
A list containing:
\item{sizes}{an integer vector of fragment lengths for all valid read pairs in the library}
\item{diagnostics}{an integer vector containing the total number of reads, the number of mapped reads, number of mapped singleton reads, pairs with exactly one unmapped read,
	number of improperly orientated read pairs and interchromosomal pairs}
}

\details{
This function computes a number of diagnostics for paired-end data in a supplied BAM file.
The aims is to provide an indication of the quality of library preparation and sequencing.

Firstly, a read is only considered to be mapped if it is not removed by \code{dedup}, \code{minq}, \code{restrict} or \code{discard} in \code{\link{readParam}}.
Otherwise, the alignment is not considered to be reliable.
Any read pair with exactly one unmapped read is discarded, and the number of read pairs lost in this manner is recorded.
Obviously, read pairs with both reads unmapped will be ignored completely, as will any unpaired reads in the BAM file.
Secondary and supplementary alignments are ignored completely and do not contribute to the total - see \code{\link{readParam}} for details.

Of the mapped pairs, the valid (i.e., proper) read pairs are identified.
This involves several criteria:
\itemize{
\item Read pairs must be intrachromosomal.
If the reads are on different chromosomes, the read pair will be recorded as being interchromosomal.
\item The two reads in the pair must lie on opposite strands.
Otherwise, the read pair will be considered as being improperly orientated.
\item The 5' end of the forward read must not map to a higher genomic coordinate than the 5' end of the reverse read.
Otherwise, the read pair will be considered as being improperly orientated.
}
Note that the 3' end of one read is allowed to overrun the 5' end of the other.
This avoids being too stringent in the presence of sequencing errors, untrimmed adaptors, etc. at the start or ends of reads.

Each valid read pair corresponds to a DNA fragment where both ends are sequenced. 
The size of the fragment can be determined by calculating the distance between the 5' ends of the mapped reads. 
The distribution of sizes is useful for assessing the quality of the library preparation, along with all of the recorded diagnostics. 
Note that any \code{max.frag} specification in \code{param} will be ignored; sizes for all valid pairs will be returned.
}

% Note; length(out$sizes)*2 == out$diagnostics[2] - sum(out$diagnostics[3:4]) - sum(out$diagnostics[5:6])*2. Because, read pairs with one mapped = 1 read in total.

\author{Aaron Lun}

\seealso{
\code{\link{readParam}}
}

\examples{
bamFile <- system.file("exdata", "pet.bam", package="csaw")
out <- getPESizes(bamFile, param=readParam(pe="both"))
out <- getPESizes(bamFile, param=readParam(pe="both", restrict="chrA"))
out <- getPESizes(bamFile, param=readParam(pe="both", discard=GRanges("chrA", IRanges(1, 50))))
}

\keyword{diagnostics}
