\name{profileSites}
\alias{profileSites}

\title{Profile binding sites}
\description{Get the coverage profile around potential binding sites.}

\usage{
profileSites(bam.files, regions, range=5000, ext=100, average=TRUE,
    weight=1, param=readParam(), strand=c("ignore", "use", "match"))
}

\arguments{
	\item{bam.files}{a character vector containing paths to one or more BAM files}
	\item{regions}{a \code{GRanges} object over which profiles are to be aggregated}
	\item{range}{an integer scalar specifying the range over which the profile will be collected}
	\item{ext}{an integer scalar specifying the average fragment length for single-end data}
	\item{average}{a logical scalar specifying whether the profiles should be averaged across regions}
	\item{weight}{a numeric vector indicating the relative weight to be assigned to the profile for each region}
	\item{param}{a \code{readParam} object containing read extraction parameters, or a list of such objects (one for each BAM file)}
	\item{strand}{a string indicating how stranded \code{regions} should be handled}
}

\details{
This function aggregates the coverage profile around the specified regions.
The shape of this profile can guide an intelligent choice of the window size in \code{\link{windowCounts}}, or to determine if region expansion is necessary in \code{\link{regionCounts}}.
For the former, restricting the \code{regions} to locally maximal windows with \code{\link{findMaxima}} is recommended prior to use of \code{profileSites}.
The function can be also used to examine average coverage around known features of interest, like genes.

The profile records the number of fragments overlapping each base within \code{range} of the start of all \code{regions}.
Single-end reads are directionally extended to \code{ext} to impute the fragment (see \code{\link{windowCounts}} for more details).
For paired-end reads, the interval between each pair is used as the fragment.
If multiple \code{bam.files} are specified, reads are pooled across files for counting into each profile.

Direct aggregation will favor high-abundance regions as these have higher counts.
If this is undesirable, high-abundance regions can be downweighted using the \code{weight} argument.
For example, this can be set to the inverse of the sum of counts across all libraries for each region in \code{regions}.
This will ensure that each region contributes equally to the final profile.

Aggregation can be turned off by setting \code{average=FALSE}. 
In such cases, a separate profile will be returned for each region, instead of the profiles being averaged across all regions.
This may be useful, e.g., for constructing heatmaps of enrichment across many regions.
Note that \code{weight} has no effect when aggregation is turned off.
}

\section{Comments on strand specificity}{
If \code{strand="use"}, the function is called separately on the reverse-stranded regions.
The profile for these regions is computed such that the left side of the profile corresponds to the upstream flank on the reverse strand (i.e., the profile is flipped).
The center of the profile corresponds to the 5' end of the region on the reverse strand.
This may be useful for features where strandedness is important, e.g., TSS's.
Otherwise, if \code{strand="ignore"}, no special treatment is given to reverse-stranded features.

By default, the strandedness of the region has no effect on read extraction.
If \code{strand="match"}, the profile for reverse-strand regions is made with reverse-strand reads only (this profile is also flipped, as described for \code{strand="use"}).
Similarly, only forward-strand reads are used for forward- or unstranded regions.
Note that \code{param$forward} must be set to \code{NULL} for this to work.
}

\value{
If \code{average=TRUE}, a numeric vector of average coverages for each base position within \code{range} is returned, where the average is taken over all \code{regions}.
The vector is named according to the relative position of each base to the start of the region.
If \code{weight} is set as described above for each region, then the vector will represent average relative coverages, i.e., relative to the number of fragments counted in the region itself.

If \code{average=FALSE}, an integer matrix of coverage values is returned. 
Each row of the matrix corresponds to an entry in \code{regions}, while each column corresponds to a base position with \code{range}.
Column names are set o the relative position of each base to the start of each region.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{findMaxima}},
\code{\link{windowCounts}},
\code{\link{wwhm}}
}

\examples{
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
data <- windowCounts(bamFile, filter=1)
rwsms <- rowSums(assay(data))
maxed <- findMaxima(rowRanges(data), range=100, metric=rwsms)
	
x <- profileSites(bamFile, rowRanges(data)[maxed], range=200)
plot(as.integer(names(x)), x)

x <- profileSites(bamFile, rowRanges(data)[maxed], range=500)
plot(as.integer(names(x)), x)

x <- profileSites(bamFile, rowRanges(data)[maxed], range=500, weight=1/rwsms)
plot(as.integer(names(x)), x)

# Introducing some strandedness.
regs <- rowRanges(data)[maxed]
strand(regs) <- sample(c("-", "+", "*"), sum(maxed), replace=TRUE)
x <- profileSites(bamFile, regs, range=500)
plot(as.integer(names(x)), x)
x2 <- profileSites(bamFile, regs, range=500, strand="use")
points(as.integer(names(x2)), x2, col="red")
x3 <- profileSites(bamFile, regs, range=500, strand="match",
    param=readParam(forward=NULL))
points(as.integer(names(x3)), x3, col="blue")

# Returning separate profiles.
y <- profileSites(bamFile, rowRanges(data)[maxed], range=500, average=FALSE)
dim(y)
}

\keyword{diagnostics}
