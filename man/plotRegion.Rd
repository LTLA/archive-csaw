\name{plotRegion}
\alias{plotRegion}

\title{Plot read coverage within a region}

\description{Plot the depth on each strand for all reads from a BAM file in a specified genomic region}

\usage{
plotRegion(cur.region, bam.file, dedup=FALSE, minq=0, pet=c("none", "both", "first", "second"), max.frag=500,
	max.depth=NULL, fcol = "blue", rcol="red", xlab="Genomic position (bp)", ylab="Read depth", ...)
}

\arguments{
	\item{cur.region}{a \code{GRanges} object of length 1 describing the region of interest}
	\item{bam.file}{a character string containing the path to a sorted and indexed BAM file}
	\item{dedup}{a logical scalar specifying whether duplicates should be removed}
    \item{minq}{an integer scalar defining the minimum threshold for the mapping quality score}	
    \item{pet}{a character string specifying whether the data involves paired end tags}
	\item{max.frag}{an integer scalar, specifying the maximum allowable fragment size for PET data}
	\item{max.depth}{a numeric value specifying the maximum height of the plot}
	\item{fcol}{a R colour object for forward reads}
	\item{rcol}{a R colour object for reverse reads}
	\item{xlab}{a character string for the x-axis label}
	\item{ylab}{a character string for the y-axis label}
	\item{...}{other arguments to be passed to \code{plot}}
}	

\details{
This function plots the read coverage on both strands for a genomic region. Reads mapping to the forward strand
are in blue and are represented as positive coverage. Reads mapping to the reverse strand are in red and
are represented as negative coverage.

If \code{pet="both"}, only valid read pairs are plotted. Coverage is plotted by considering the
fragment represented by each read pair and splitting that fragment in two. The 5' half is plotted as a forward
read and the 3' half is plotted as the reverse read.

The maximum depth should be scaled by the relative library size. For example, smaller libraries should have a
smaller \code{max.depth} to account for the reduced total number of reads. If \code{max.depth=NULL}, the ideal 
height of the plot is estimated from the maximum base coverage within the specified region.

As always, the \code{dedup} argument only works if duplicate reads have been marked with the appropriate flag in
the BAM file.
}

\value{
A plot is produced on the current graphics device. A \code{NULL} value is invisibly returned.
}

\author{Aaron Lun}

\examples{
bamFiles<-system.file("exdata", "chip", c("rep1_sort.bam", "rep2_sort.bam"), package="csaw")

# Normalization factors are set to 1, so library sizes are unaffected.
lib.sizes<-c(countBam(bamFiles[1])$records, countBam(bamFiles[2])$records)
mlib<-mean(lib.sizes)

# Plotting to an adjusted max.depth for (effective) library sizes. 
constant<-5
par(mfrow=c(1,2))
plotRegion(GRanges("chrA", IRanges(100, 500)), bamFiles[1], max.depth=constant*lib.sizes[1]/mlib)
plotRegion(GRanges("chrA", IRanges(100, 500)), bamFiles[2], max.depth=constant*lib.sizes[2]/mlib)

# Considering paired-read data.
bamFile <-system.file("exdata", "pet", c("pet_sort.bam"), package="csaw")
plotRegion(cur.region=GRanges("chrA", IRanges(0, 200)), bam.file=bamFile, pet="both")
plotRegion(cur.region=GRanges("chrA", IRanges(0, 200)), bam.file=bamFile, pet="both", max.frag=50)
}

\references{
Carlson M and Aboyoun P (2013). GenomicRanges Use Cases. \emph{Bioconductor} (version 2.11)
}
