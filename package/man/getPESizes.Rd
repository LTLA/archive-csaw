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
This function assembles a number of paired-end diagnostics.
For starters, a read is only mapped if it is not removed by \code{dedup}, \code{minq}, \code{restrict} or \code{discard} in \code{\link{readParam}}.
Otherwise, the alignment is not considered to be reliable.
Any read pair with exactly one unmapped read is discarded, and the number of read pairs lost in this manner is recorded.
Obviously, read pairs with both reads unmapped will be ignored completely.

Of the mapped pairs, the valid (i.e., proper) read pairs are identified.
These refer to intrachromosomal read pairs where the reads with the lower and higher genomic coordinates map to the forward and reverse strand, respectively.
The distance between the positions of the mapped 5' ends of the two reads must also be equal to or greater than the read lengths. 
Any intrachromosomal read pair that fails these criteria will be considered as improperly oriented.
If the reads are on different chromosomes, the read pair will be recorded as being interchromosomal.

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
