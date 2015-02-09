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
Users should be aware that the downstream range-processing functions are generally not strand-aware, e.g., \code{\link{mergeWindows}}.
Any strandedness of the ranges will be ignored in these functions.
If strand-specific processing is desired, users must manually subset the input by strand, and then reassemble the function outputs - see the examples below.
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
is.forward <- as.logical(strand(out) == "+")
fclust <- mergeWindows(rowData(out[is.forward]), tol=100)
rclust <- mergeWindows(rowData(out[!is.forward]), tol=100)

strand(fclust$region) <- "+"
strand(rclust$region) <- "-"
cluster.ids <- integer(nrow(out))
cluster.ids[is.forward] <- fclust$id
cluster.ids[!is.forward] <- rclust$id + length(fclust$region)
clustered <- list(id=cluster.ids, region=c(fclust$region, rclust$region))

# Running findMaxima with strand-specificity.
rwsms <- rowSums(assay(out))
f.ismax <- findMaxima(rowData(out[is.forward]), range=100, metric=rwsms[is.forward])
r.ismax <- findMaxima(rowData(out[!is.forward]), range=100, metric=rwsms[!is.forward])
ismax <- logical(nrow(out))
ismax[is.forward] <- f.ismax
ismax[!is.forward] <- r.ismax
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}}
}

\author{
Aaron Lun
}
