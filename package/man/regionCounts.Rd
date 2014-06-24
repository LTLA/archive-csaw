\name{regionCounts}
\alias{regionCounts}

\title{Count reads overlapping each region}

\description{Count the number of extended reads overlapping pre-specified regions}

\usage{
regionCounts(bam.files, regions, ext=100, 
    pet=c("none", "both", "first", "second"), max.frag=500, 
    rescue.pairs=FALSE, dedup=FALSE, minq=NA, restrict=NULL, 
    discard=NULL)
}

\arguments{
	\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
	\item{regions}{a \code{GRanges} object containing the regions over which reads are to be counted}
	\item{ext}{an integer scalar describing the average length of the sequenced fragment}
	\item{pet}{a character string specifying whether the data involves paired end tags}
	\item{max.frag}{an integer scalar, specifying the maximum allowable fragment size for PET data}
	\item{rescue.pairs}{a logical scalar indicating whether invalid read pairs should be counted}
	\item{dedup}{a boolean specifying whether marked duplicates should be removed}
	\item{minq}{an integer scalar defining the minimum threshold for the mapping quality score, none by default}
    \item{restrict}{a character vector specifying which chromosomes should be used for read counting}
	\item{discard}{a \code{GRanges} object specifying the regions in which reads are to be ignored}
}

\value{
A list is returned containing:
	\item{counts}{an integer matrix containing counts for each library at each region}
    \item{totals}{an integer vector containing the total number of reads in each library}
}

\details{
This function simply provides a wrapper around \code{\link{countOverlaps}}
for read counting into specified regions. It is provided so as to allow
for counting with awareness of the other parameters, e.g., \code{ext},
\code{pet}, in case users wish to coordinate these counts with counts from
\code{\link{windowCounts}}. Checking that the output \code{totals} are
the same is strongly recommended.
}

\author{Aaron Lun}

\seealso{
\code{\link{countOverlaps}},
\code{\link{windowCounts}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
incoming <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(200, 1000, 700, 1500)))
regionCounts(bamFiles, regions=incoming)
regionCounts(bamFiles, regions=incoming, restrict="chrB")

# Loading PET data.
bamFile <- system.file("exdata", "pet.bam", package="csaw")
regionCounts(bamFile, pet="both", regions=incoming)
regionCounts(bamFile, max.frag=100, pet="first", 
    regions=incoming, restrict="chrA")
regionCounts(bamFile, max.frag=100, pet="both", 
    regions=incoming, restrict="chrA", rescue.pairs=TRUE)
}

\keyword{counting}
