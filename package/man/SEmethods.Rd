\name{SEmethods}
\alias{normalize}
\alias{normalize,RangedSummarizedExperiment-method}
\alias{normOffsets,RangedSummarizedExperiment-method}
\alias{asDGEList,RangedSummarizedExperiment-method}
\alias{asDGEList}

\title{Statistical wrappers for RangedSummarizedExperiment objects}

\description{Convenience wrappers for statistical routines operating on RangedSummarizedExperiment objects.}

\usage{
\S4method{normOffsets}{RangedSummarizedExperiment}(object, lib.sizes, ...)
\S4method{normalize}{RangedSummarizedExperiment}(object, lib.sizes, ...) # deprecated, use normOffsets
\S4method{asDGEList}{RangedSummarizedExperiment}(object, lib.sizes, ...)
}

\arguments{
  \item{object}{a \code{RangedSummarizedExperiment} object, like that produced by \code{\link{windowCounts}}}
  \item{lib.sizes}{an (optional) integer vector of library sizes}
  \item{...}{other arguments to be passed to the function being wrapped}
}

\details{
Counts are extracted from the matrix corresponding to the first assay in the \code{RangedSummarizedExperiment} object.
If not specified in \code{lib.sizes}, library sizes are taken from the \code{totals} field in the column data of \code{object}.
Warnings will be generated if this field is not present.

In the \code{normOffsets} method, the extracted counts and library sizes are supplied to \code{\link{normOffsets,matrix-method}}.
Similarly, the \code{asDGEList} method wraps the \code{\link{DGEList}} constructor.
In both cases, any arguments in \code{...} are also passed to the wrapped functions.
}

\value{
For \code{normOffsets}, either a numeric matrix or vector is returned; see \code{\link{normOffsets,matrix-method}}.

For \code{asDGEList}, a \code{DGEList} object is returned.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)
normOffsets(data)
normOffsets(data, lib.sizes=c(10, 100))
head(normOffsets(data, type="loess"))

asDGEList(data)
asDGEList(data, lib.sizes=c(10, 100))
asDGEList(data, norm.factors=c(1.11, 2.23), group=c("a", "b"))
}

\seealso{
\code{\link{normOffsets}},
\code{\link{DGEList}},
\code{\link{windowCounts}}
}

\keyword{normalization}
