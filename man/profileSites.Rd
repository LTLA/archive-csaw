\name{profileSites}
\alias{profileSites}

\title{Profile binding sites}
\description{Get the coverage profile around potential binding sites.}

\usage{
profileSites(bam.files, regions, param=readParam(), range=5000, ext=100, 
    average=TRUE, normalize="none", strand=c("ignore", "use", "match"))
}

\arguments{
	\item{bam.files}{a character vector containing paths to one or more BAM files}
	\item{regions}{a \code{GRanges} object over which profiles are to be aggregated}
	\item{param}{a \code{readParam} object containing read extraction parameters}
	\item{range}{an integer scalar specifying the range over which the profile will be collected}
	\item{ext}{an integer scalar or list specifying the average fragment length for single-end data}
	\item{average}{a logical scalar specifying whether the profiles should be averaged across regions}
    \item{normalize}{a string specifying how normalization of each region's profile should be performed prior to averaging}
	\item{strand}{a string indicating how stranded \code{regions} should be handled}
}

\details{
This function computes the average coverage profile around the specified regions.
Specifically, the profile is constructed by counting the number of fragments overlapping each base in the interval flanking each entry of \code{regions}. 
The interval for each entry is centred at its start location (base zero) and spans the flanking \code{range} on either side.

Single-end reads are directionally extended to \code{ext} to impute the fragment (see \code{\link{windowCounts}} for more details).
For paired-end reads, the interval between each pair is used as the fragment.
If multiple \code{bam.files} are specified, reads are pooled across files for counting into each profile.

By default, an average of the coverage profiles across all \code{regions} is returned.
Normalization of each region's profile is performed on by setting \code{normalize} to:
\describe{
\item{\code{none}:}{No normalization is performed, i.e., counts per base are directly averaged across all regions.
Thus, the shape of the average profile is largely determined by high-abundance regions.
}
\item{\code{total}:}{The profile for each region is divided by the sum of coverages across all bases in the interval.
This effectively normalizes for the total number of reads in each region.
}
\item{\code{max}:}{The profile for each region is divided by its maximum value.
This ensures that the maximum height of each region is the same.
} 
}
If \code{average=FALSE}, a separate profile will be returned for each region instead.
This may be useful, e.g., for constructing heatmaps of enrichment across many regions.

The profile can be used to examine average coverage around known features of interest, like genes or transcription start sites.
Its shape can guide the choice of the window size in \code{\link{windowCounts}}, or to determine if larger regions should be used in \code{\link{regionCounts}}.
For the former, restricting the \code{regions} to locally maximal windows with \code{\link{findMaxima}} is recommended to capture the profile of binding events.
}

\section{Comments on strand specificity}{
By default, the strandedness of the regions are ignored with \code{strand="ignore"}.
If \code{strand="use"}, the behaviour of this function will differ between forward- and reverse-stranded entries in \code{regions}.
\itemize{
\item Forward-stranded or unstranded regions are processed as previously described above.
Base zero corresponds to the start of the region, negative distances correspond to the 5' flanking region, and positive distances correspond to the 3' flanking region.
\item Reverse-stranded regions are flipped, i.e., base zero corresponds to the \emph{end} of the region.
Negative distances correspond to the 5' flanking region on the reverse strand, while positive distances correspond to the 3' flanking region on this strand.
}
This ensures that the center of the profile always corresponds to the 5' end of the region, with upstream regions on the left and downstream regions on the right.
This may be useful for features where strandedness is important, e.g., TSS's.

By default, the strandedness of the region has no effect on the choice of reads that are used.
If \code{strand="match"}, the profile for reverse-strand regions is constructed from reverse-strand reads only.
Similarly, only forward-strand reads are used for forward- or unstranded regions.
Note that \code{param$forward} must be set to \code{NULL} for this to work.
Flipping of reverse-stranded profiles is also performed in this setting, as described for \code{strand="use"}.
}

\value{
If \code{average=TRUE}, a numeric vector of average coverages for each base position within \code{range} is returned, where the average is taken over all \code{regions}.
The vector is named according to the relative position of each base to the start of the region.
The interpretation of the coverages will depend on the value of \code{normalize}.

If \code{average=FALSE}, an integer matrix of coverage values is returned. 
Each row of the matrix corresponds to an entry in \code{regions}, while each column corresponds to a base position with \code{range}.
Column names are set to the relative position of each base to the start of each region.
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

# Running profileSites .
x <- profileSites(bamFile, rowRanges(data)[maxed], range=200)
plot(as.integer(names(x)), x)

x <- profileSites(bamFile, rowRanges(data)[maxed], range=500)
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
