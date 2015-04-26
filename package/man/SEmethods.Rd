\name{SEmethods}
\alias{normalize,SummarizedExperiment-method}
\alias{normalize}
\alias{asDGEList,SummarizedExperiment-method}
\alias{asDGEList}

\title{Statistical wrappers for SummarizedExperiment objects}

\description{Convenience wrappers for statistical routines operating on SummarizedExperiment objects}

\usage{
normalize(object, ...)
asDGEList(object, ...)
}

\arguments{
  \item{object}{a \code{SummarizedExperiment} object, like that produced by \code{\link{windowCounts}}}
  \item{...}{other arguments to be passed to the function being wrapped}
}

\details{
Counts are extracted from the matrix corresponding to the first assay in the \code{SummarizedExperiment} object.
The total library size is taken from the \code{totals} field in the column data of \code{object}.
Warnings will be generated if this field is not present.

In the \code{normalize} method, the extracted counts and library sizes are supplied to \code{\link{normalizeCounts}}, along wth arguments in \code{...}.
Similarly, the \code{asDGEList} method wraps the \code{\link{DGEList}} constructor.
Library sizes can be manually specified in \code{...} by setting \code{lib.sizes} for \code{normalize}, or \code{lib.size} for \code{asDGEList}.
}

\value{
For \code{normalize}, either a numeric matrix or vector is returned; see \code{\link{normalizeCounts}}.

For \code{asDGEList}, a \code{DGEList} object is returned.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)
normalize(data)
normalize(data, lib.sizes=c(10, 100))
head(normalize(data, type="loess"))

asDGEList(data)
asDGEList(data, lib.size=c(10, 100))
asDGEList(data, norm.factors=c(1.11, 2.23), group=c("a", "b"))
}

\seealso{
\code{\link{normalizeCounts}},
\code{\link{DGEList}},
\code{\link{windowCounts}}
}

\keyword{normalization}
