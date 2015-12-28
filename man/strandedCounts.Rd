\name{strandedCounts}
\alias{strandedCounts}

\title{Get strand-specific counts}
\description{Obtain strand-specific counts for each genomic window or region.}

\usage{
strandedCounts(bam.files, param=readParam(forward=NULL), regions=NULL, ...)
}

\arguments{
\item{bam.files}{a character vector containing paths to sorted and indexed BAM files}
\item{param}{a \code{readParam} object containing read extraction parameters, or a list of such objects (one for each BAM file), where the \code{forward} slot must be set to \code{NULL}}
\item{regions}{a \code{GRanges} object specifying the regions over which reads are to be counted}
\item{...}{other arguments to be passed to \code{\link{windowCounts}} or \code{\link{regionCounts}}}
}

\value{
A \code{RangedSummarizedExperiment} object containing strand-specific counts for genomic regions.
}

\details{
Some applications require strand-specific counts for each genomic region.
This function calls \code{\link{windowCounts}} after setting \code{param$forward} to \code{TRUE} and \code{FALSE}.
Any existing value of \code{param$forward} is ignored.
If \code{regions} is specified, \code{\link{regionCounts}} is used instead of \code{\link{windowCounts}}.

The function then concatenates the two \code{RangedSummarizedExperiment} objects (one from each strand).
The total numbers of reads are added together to form the new \code{totals}.
However, the total numbers of reads for each strand are also stored for future reference.
Count loading parameters are also stored in the \code{metadata}.

Each row in the concatenated object corresponds to a stranded genomic region, where the strand of the region indicates the strand of the reads that were counted in that row.
Note that there may not be two rows for each genomic region.
This is because any empty rows, or those with counts below \code{filter}, will be removed within each call to \code{\link{windowCounts}}.
}

\section{Warnings}{
Users should be aware that many of the downstream range-processing functions are not strand-aware by default, e.g., \code{\link{mergeWindows}}.
Any strandedness of the ranges will be ignored in these functions.
If strand-specific processing is desired, users must manually set \code{ignore.strand=FALSE}.

The input \code{param$forward} should be set to \code{NULL}, as a safety measure.
This is because the returned object is a composite of two separate calls to the relevant counting function.
If the same \code{param} object is supplied to other functions, an error will be thrown if \code{param$forward} is \code{NULL}.
This serves to remind users that such functions should instead be called twice, i.e., separately for each strand after setting \code{param$forward} to \code{TRUE} or \code{FALSE}.
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
xparam <- readParam(forward=NULL)
out <- strandedCounts(bamFiles, filter=1, param=xparam)
out

strandedCounts(bamFiles, filter=1, width=100, param=xparam)
strandedCounts(bamFiles, filter=1, param=reform(xparam, minq=20))

incoming <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(200, 1000, 700, 1500)))
strandedCounts(bamFiles, regions=incoming, param=xparam)
strandedCounts(bamFiles, regions=incoming, param=reform(xparam, dedup=TRUE))

# Throws an error, as the same reads are not involved.
try(windowCounts(bamFiles, filter=1, width=100, param=xparam))

# Library sizes should be the same as that without strand-specificity.
colData(out)
out.ref <- windowCounts(bamFiles, param=reform(xparam, forward=NA))
stopifnot(identical(out.ref$totals, out$totals))

# Running assorted functions on strandedCounts output.
mergeWindows(rowRanges(out), tol=100)
mergeWindows(rowRanges(out), tol=100, ignore.strand=FALSE)

rwsms <- rowSums(assay(out))
summary(findMaxima(rowRanges(out), range=100, metric=rwsms))
summary(findMaxima(rowRanges(out), range=100, metric=rwsms, ignore.strand=FALSE))
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}}
}

\author{
Aaron Lun
}

\keyword{counting}
