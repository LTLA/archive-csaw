\name{Parameter list methods}
\alias{param,SummarizedExperiment-method}
\alias{param}
\alias{reformList}

\title{Get or modify readParam lists}
\description{Extract and modify lists of library-specific readParam objects}

\usage{
param(object)
reformList(paramlist, ...)
}

\arguments{
  \item{object}{a \code{SummarizedExperiment} object, like that produced by \code{\link{windowCounts}}}
  \item{paramlist}{a list of \code{readParam} objects}
  \item{...}{other arguments, to be passed to \code{\link{reform}}}
}

\details{
The \code{param} function extracts a \code{readParam} object, or a list of such objects, from a \code{SummarizedExperiment} object.
The former is returned if a common \code{readParam} was used for counting in \code{\link{windowCounts}} or \code{\link{regionCounts}}.
Otherwise, the latter is returned if a different \code{readParam} was used for each library.

The \code{reformList} function sets the parameters specified in \code{...} across all elements of \code{paramlist}.
This standardizes all of the internal \code{readParam} objects, e.g., for the \code{restrict} parameter.
If \code{paramlist} is a single \code{readParam} object, then \code{reformList} is equivalent to \code{\link{reform}}.
}

\value{
For both functions, a common \code{readParam} object or a list of library-specific \code{readParam} objects.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)
param(data)

data <- windowCounts(bamFiles, width=100, filter=1, param=list(readParam(), readParam(minq=100)))
param(data)
param(data[,1])
param(data[,2])

plist <- param(data)
reformList(plist, restrict="chrA")
reformList(plist, minq=10)
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}},
\code{\link{reform}}
}

