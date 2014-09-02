\name{normalize}
\alias{normalize,SummarizedExperiment-method}
\alias{normalize}

\title{Normalize counts between ChIP-seq libraries}

\description{Calculate normalization factors or offsets for multiple ChIP-seq libraries.}

\usage{
normalize(object, ...)
}

\arguments{
  \item{object}{a \code{SummarizedExperiment} object, like that produced by \code{\link{windowCounts}}}
  \item{...}{other arguments to be passed to \code{\link{normalizeCounts}}}
}

\details{
This method provides a convenience wrapper around \code{\link{normalizeCounts}} for \code{SummarizedExperiment} objects.
If \code{counts} is a \code{SummarizedExperiment} object, the counts are extracted using the matrix corresponding to the first assay.
The total library size is taken from the \code{totals} entry in the column data.
Appropriate warnings will be supplied if this entry is not available.
}

\value{
Either a numeric matrix or vector; see \code{\link{normalizeCounts}}.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)
normalize(data)
head(normalize(data, type="loess"))
}

\seealso{
\code{\link{normalizeCounts}},
\code{\link{windowCounts}}
}

\keyword{normalization}
