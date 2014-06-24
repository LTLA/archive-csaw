\name{getPETSizes}
\alias{getPETSizes}

\title{Compute fragment lengths for paired-end tags}
 
\description{Compute the length of the sequenced fragment for each read pair in paired-end tag (PET) data.}

\usage{
getPETSizes(bam.file, dedup=FALSE, minq=NA, restrict=NULL, 
    discard=NULL)
}

\arguments{
	\item{bam.file}{a character string containing the file path to a sorted and indexed BAM file}
	\item{dedup}{a boolean scalar specifying whether or not marked duplicates should be removed}
    \item{minq}{an integer scalar defining the minimum threshold for the mapping quality score when computing sizes, none by default}
    \item{restrict}{a character vector specifying chromosomes from which fragment sizes should be extracted}
    \item{discard}{a \code{GRanges} object specifying the regions in which reads are to be ignored}
}

\value{
A list containing:
\item{sizes}{an integer vector of fragment lengths for all valid read pairs in the library}
\item{diagnostics}{an integer vector containing the total number of reads, number of singleton reads, number of improperly orientated read pairs,
	pairs with one unmapped read and interchromosomal pairs}
}

\details{
Each read pair corresponds to a DNA fragment where both ends are sequenced. The size of the fragment can be determined by the distance between the 5' ends 
of the mapped reads. The distribution of sizes is useful for assessing the quality of the library preparation. Some diagnostics are also returned to 
provide some indication of the quality of the paired-end sequencing and alignment.
		
Only the sizes of valid read pairs are considered. The read with the lower genomic coordinate must map to the forward strand and the read with the higher 
genomic coordinate must map to the reverse strand. The distance between the positions of the mapped 5' ends of each read must also be equal to or greater 
than the read lengths. Otherwise, the read pair will be considered as improperly orientated.

The combined effects of \code{dedup}, \code{minq}, \code{restrict} and \code{discard} will define the total number of reads that are considered mapped.
Read pairs must pass all criteria and be successfully mapped on two separate chromosomes in order to be considered as interchromosomal. Otherwise, 
the read pair will be recorded as having one unmapped mate as the alignment is not reliable if it fails any of the criteria.

Note that any duplicate read must be marked in the bit field of the BAM file using a tool like Picard if it is to be removed with \code{dedup=TRUE}.
Also note that the read pair information in each BAM file \emph{must} be synchronized.
}

\author{Aaron Lun}

\seealso{
\code{\link{windowCounts}}
}

\examples{
bamFile <- system.file("exdata", "pet.bam", package="csaw")
out <- getPETSizes(bamFile)
out <- getPETSizes(bamFile, restrict="chrA")
out <- getPETSizes(bamFile, discard=GRanges("chrA", IRanges(1, 50)))
}

\keyword{diagnostics}
