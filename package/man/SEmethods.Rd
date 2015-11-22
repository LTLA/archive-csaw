\name{SEmethods}
\alias{normalize}
\alias{normalize,SummarizedExperiment0-method}
\alias{normOffsets,SummarizedExperiment0-method}
\alias{asDGEList,SummarizedExperiment0-method}
\alias{asDGEList}

\title{Statistical wrappers for SummarizedExperiment0 objects}

\description{Convenience wrappers for statistical routines operating on SummarizedExperiment0 objects.}

\usage{
\S4method{normOffsets}{SummarizedExperiment0}(object, lib.sizes, assay=1, ...)
\S4method{normalize}{SummarizedExperiment0}(object, lib.sizes, type="scaling", assay=1, ...) 
\S4method{asDGEList}{SummarizedExperiment0}(object, lib.sizes, norm.factors, ...)
}

\arguments{
  \item{object}{a SummarizedExperiment0 object or its derived classes, like that produced by \code{\link{windowCounts}}}
  \item{lib.sizes}{an (optional) integer vector of library sizes}
  \item{assay}{a numeric scalar indicating which assay in \code{object} contains the count matrix}
  \item{type}{a character scaling indicating the type of normalization, see \code{\link{normOffsets}}}
  \item{norm.factors}{an (optional) numeric vector of normalization factors}
  \item{...}{other arguments to be passed to the function being wrapped}
}

\details{
For \code{normalize} and \code{normOffset}, counts are extracted from the matrix corresponding to the first assay in the SummarizedExperiment0 object 
    (different assays can be specified with the \code{assay} argument).
For \code{asDGEList}, counts are extracted from the assay named \code{"counts"}, with the first assay used if no such assay exists with that name.
If not specified in \code{lib.sizes}, library sizes are taken from the \code{totals} field in the column data of \code{object}.
Warnings will be generated if this field is not present.

In the \code{normOffsets} and \code{normalize} methods, the extracted counts and library sizes are supplied to \code{\link{normOffsets,matrix-method}}.
Similarly, the \code{asDGEList} method wraps the \code{\link{DGEList}} constructor.
In both cases, any arguments in \code{...} are also passed to the wrapped functions.

If \code{type="scaling"} in \code{normalize}, a vector of normalization factors is stored in the \code{norm.factors} field of the \code{colData} in the output object.
If \code{type="loess"}, an offset matrix is instead stored as a named entry in the \code{assays} of the output object.

If \code{norm.factors} is not specified, \code{asDGEList} will attempt to extract normalization factors from \code{object$norm.factors}.
If this is not available, factors will be set to the default (all unity).
If \code{assays(object)$offset} is present, this will be assigned to the \code{offset} field of the output \code{DGEList} object.
}

\value{
For \code{normOffsets}, either a numeric matrix or vector is returned; see \code{\link{normOffsets,matrix-method}}.

For \code{normalize}, an object of the same class as \code{object} is returned with normalization data stored in various fields.

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

norm.data <- normalize(data)
norm.data$norm.factors
asDGEList(norm.data)$samples$norm.factors
asDGEList(norm.data, norm.factors=c(1,2))$samples$norm.factors # gets priority
norm.data <- normalize(data, type="loess")
head(asDGEList(norm.data)$offset)
}

\seealso{
\code{\link{normOffsets}},
\code{\link{DGEList}},
\code{\link{windowCounts}}
}

\keyword{normalization}
