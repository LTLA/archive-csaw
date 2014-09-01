\name{getPETSizes}
\alias{getPETSizes}

\title{Compute fragment lengths for paired-end tags}
 
\description{Compute the length of the sequenced fragment for each read pair in paired-end tag (PET) data.}

\usage{
getPETSizes(bam.file, param=readParam(pet="both"))
}

\arguments{
	\item{bam.file}{a character string containing the file path to a sorted and indexed BAM file}
    \item{param}{a \code{readParam} object containing read extraction parameters}
}

\value{
A list containing:
\item{sizes}{an integer vector of fragment lengths for all valid read pairs in the library}
\item{diagnostics}{an integer vector containing the total number of reads, number of singleton reads, number of improperly orientated read pairs,
	pairs with one unmapped read and interchromosomal pairs}
}

\details{
Each read pair corresponds to a DNA fragment where both ends are sequenced. 
The size of the fragment can be determined by the distance between the 5' ends of the mapped reads. 
The distribution of sizes is useful for assessing the quality of the library preparation. 
Some diagnostics are also returned to provide some indication of the quality of the paired-end sequencing and alignment.
		
Only the sizes of valid read pairs are considered.
In valid/proper pairs, the read with the lower genomic coordinate must map to the forward strand and the read with the higher genomic coordinate must map to the reverse strand. 
The distance between the positions of the mapped 5' ends of each read must also be equal to or greater than the read lengths. 
Otherwise, the read pair will be considered as improperly orientated.

The combined effects of \code{dedup}, \code{minq}, \code{restrict} and \code{discard} in \code{\link{readParam}} will define the total number of reads that are considered mapped.  
Read pairs must pass all criteria and be successfully mapped on two separate chromosomes in order to be considered as interchromosomal. 
Otherwise, the read pair will be recorded as having one unmapped mate. 
An alignment is not considered to be reliable if it fails any of the criteria.
}

\author{Aaron Lun}

\seealso{
\code{\link{readParam}}
}

\examples{
bamFile <- system.file("exdata", "pet.bam", package="csaw")
out <- getPETSizes(bamFile, param=readParam(pet="both"))
out <- getPETSizes(bamFile, param=readParam(pet="both", restrict="chrA"))
out <- getPETSizes(bamFile, param=readParam(pet="both", discard=GRanges("chrA", IRanges(1, 50))))
}

\keyword{diagnostics}
