\name{getPETSizes}
\alias{getPETSizes}

\title{Compute fragment lengths for paired-end tags}
 
\description{Compute the length of the sequenced fragment for each read pair in paired-end tag (PET) data.}

\usage{
getPETSizes(bam.file, dedup=FALSE, minq=0, restrict=NULL)
}

\arguments{
	\item{bam.file}{a character string containing the file path to a sorted and indexed BAM file}
	\item{dedup}{a boolean scalar specifying whether or not marked duplicates should be removed}
    \item{minq}{an integer scalar defining the minimum threshold for the mapping quality score when computing sizes}
    \item{restrict}{a character vector specifying chromosomes from which fragment sizes should be extracted}
}

\value{
A list containing:
\item{sizes}{an integer vector of fragment lengths for all valid read pairs in the library}
\item{diagnostics}{an integer vector containing the total number of reads, number of singleton reads, pairs with one unmapped read and interchromosomal pairs}
}

\details{
Each read pair corresponds to a DNA fragment where both ends are sequenced. The size of the fragment can be determined by the distance between the 5' ends 
of the mapped reads. The distribution of sizes is useful for assessing the quality of the library preparation. Some diagnostics are also returned to 
provide some indication of the quality of the paired-end sequencing and alignment.
		
Only the sizes of valid read pairs are considered. The read with the lower genomic coordinate must map to the forward strand and the read with the higher 
genomic coordinate must map to the reverse strand. The distance between the positions of the mapped 5' ends of each read must also be equal to or greater 
than the read lengths. 

% Any difference in total count in \code{diagnostics} and sum of valid, interchromosomal and singleton reads represents invalid PETs which are not otherwise counted e.g. incorrect mate orientation.
		
Note that any duplicate read must be marked in the bit field of the BAM file using a tool like Picard if it is to be removed with \code{dedup=TRUE}.
Read pair information should also be properly synchronised prior to the use of this function. 

Only the mapping quality score of the read on the forward strand is used to filter read pairs if \code{minq>0}. In addition, this is only performed when 
computing fragment sizes for valid read pairs. No filtering on mapping quality is performed when computing the values returned in \code{diagnostics}.

% Doesn't really matter, as the requirement for valid pairs is usually stringent enough. It's also arguable that the MAPQ for both 
% reads in a read pair should be the same as the presence of one read should direct the uniqueness of alignment of the other read.
}

\author{Aaron Lun}

\seealso{
\code{\link{windowCounts}}
}

\examples{
bamFile<-system.file("exdata", "pet", "pet_sort.bam", package="csaw")
out<-getPETSizes(bamFile)
out<-getPETSizes(bamFile, restrict="chrA")
}

