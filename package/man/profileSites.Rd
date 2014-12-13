\name{profileSites}
\alias{profileSites}

\title{Profile binding sites}
\description{Get the coverage profile around potential binding sites.}

\usage{
profileSites(bam.files, regions, range=5000, ext=100, weight=1, final.ext=NULL, param=readParam()) 
}

\arguments{
	\item{bam.files}{a character vector containing paths to one or more BAM files}
	\item{regions}{a \code{GRanges} object over which profiles are to be aggregated}
	\item{range}{an integer scalar specifying the range over which the profile will be collected}
	\item{ext}{an integer scalar specifying the average fragment length for single-end data}
	\item{weight}{a numeric vector indicating the relative weight to be assigned to the profile for each region}
	\item{final.ext}{an integer scalar indicating the desired fragment length for extended reads}
	\item{param}{a \code{readParam} object containing read extraction parameters, or a list of such objects (one for each BAM file)}
}

\details{
This function aggregates the coverage profile around the specified regions.
The shape of this profile can guide an intelligent choice of the window size in \code{\link{windowCounts}}, or to determine if region expansion is necessary in \code{\link{regionCounts}}.
For the former, restricting the \code{regions} to locally maximal windows with \code{\link{findMaxima}} is recommended prior to use of \code{profileSites}.

More specifically, for each region, the profile consists of the number of non-overlapping fragments at each distance away from the closest edge of the region (to a maximum distance of \code{range}).
Single-end reads are directionally extended to \code{ext} to impute the fragment (with some modification according to \code{final.ext}; see \code{\link{windowCounts}} for more details).
For paired-end reads, the interval between each pair is used as the fragment.
If multiple \code{bam.files} are specified, reads are pooled across files for counting into each profile.

Direct aggregation will favour high-abundance regions as these have higher counts.
If this is undesirable, high-abundance regions can be downweighted using the \code{weight} argument.
For example, this can be set to the inverse of the sum of counts across all libraries for each region in \code{regions}.
This will ensure that each region contributes equally to the final profile.
}

\value{
A numeric vector of average coverages for each distance, where the average is taken over all \code{regions}.
If \code{weight} is set as described above for each region, then the vector will represent average relative coverages, i.e., relative to the number of fragments counted in the region itself.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{findMaxima}},
\code{\link{windowCounts}}
}

\examples{
bamFile <- system.file("exdata", "rep1.bam", package="csaw")
data <- windowCounts(bamFile, filter=1)
rwsms <- rowSums(assay(data))
maxed <- findMaxima(rowData(data), range=100, metric=rwsms)
	
x <- profileSites(bamFile, rowData(data)[maxed], range=200)
plot(x)

x <- profileSites(bamFile, rowData(data)[maxed], range=500)
plot(x)

x <- profileSites(bamFile, rowData(data)[maxed], range=500, weight=1/rwsms)
plot(x)
}

\keyword{diagnostics}
