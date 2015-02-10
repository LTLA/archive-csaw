\name{strandedCounts}
\alias{strandedCounts}

\title{Get strand-specific counts}
\description{Obtain strand-specific counts for each genomic window or region.}

\usage{
strandedCounts(bam.files, param=readParam(), regions=NULL, ...)
}

\arguments{
\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
\item{param}{a \code{readParam} object containing read extraction parameters, or a list of such objects (one for each BAM file)}
\item{regions}{a \code{GRanges} object specifying the regions over which reads are to be counted}
\item{...}{other arguments to be passed to \code{\link{windowCounts}} or \code{\link{regionCounts}}}
}

\value{
A \code{SummarizedExperiment} object containing strand-specific counts for genomic regions.
}

\details{
Some applications require strand-specific counts for each genomic region.
This function calls \code{\link{windowCounts}} after setting \code{param$forward} to \code{TRUE} and \code{FALSE}.
Any existing value of \code{param$forward} is ignored.
If \code{regions} is specified, \code{\link{regionCounts}} is used instead of \code{\link{windowCounts}}.

The function then concatenates the two \code{SummarizedExperiment} objects (one from each strand).
The total numbers of reads are added together to form the new \code{totals}.
However, the total numbers of reads for each strand are also stored for future reference.
Count loading parameters are also stored in the \code{exptData}.

Each row in the concatenated object corresponds to a stranded genomic region, where the strand of the region indicates the strand of the reads that were counted in that row.
Note that there may not be two rows for each genomic region.
This is because any empty rows, or those with counts below \code{filter}, will be removed within each call to \code{\link{windowCounts}}.
}

\section{Warning}{
Users should be aware that the downstream range-processing functions are not strand-aware by default, e.g., \code{\link{mergeWindows}}.
Any strandedness of the ranges will be ignored in these functions.
If strand-specific processing is desired, users must manually set \code{ignore.strand=FALSE}.
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
out <- strandedCounts(bamFiles, filter=1)
out
strandedCounts(bamFiles, filter=1, width=100)
strandedCounts(bamFiles, filter=1, param=readParam(minq=20))

incoming <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(200, 1000, 700, 1500)))
strandedCounts(bamFiles, regions=incoming)
strandedCounts(bamFiles, regions=incoming, param=readParam(dedup=TRUE))

# Running mergeWindows with strand-specificity.
mergeWindows(rowData(out), tol=100)
mergeWindows(rowData(out), tol=100, ignore.strand=FALSE)

rwsms <- rowSums(assay(out))
summary(findMaxima(rowData(out), range=100, metric=rwsms))
summary(findMaxima(rowData(out), range=100, metric=rwsms, ignore.strand=FALSE))
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}}
}

\author{
Aaron Lun
}
