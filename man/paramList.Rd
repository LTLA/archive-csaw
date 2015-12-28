\name{Parameter list methods}
\alias{reformList}
\alias{checkList}

\title{Get or modify readParam lists}
\description{Extract and modify lists of library-specific readParam objects}

\usage{
reformList(paramlist, ...)
checkList(paramlist)
}

\arguments{
  \item{paramlist}{a list of \code{readParam} objects}
  \item{...}{other arguments, to be passed to \code{\link{reform}}}
}

\details{
The \code{reformList} function sets the parameters specified in \code{...} across all elements of \code{paramlist}.
This standardizes all of the internal \code{readParam} objects, e.g., for the \code{restrict} parameter.
If \code{paramlist} is a single \code{readParam} object, then \code{reformList} is equivalent to \code{\link{reform}}.

The \code{checkList} function checks which parameter settings are identical between elements of \code{paramlist}.
This is useful when some parameters must be variable, while others must be the same.
}

\value{
For \code{reformList}, a list of library-specific \code{readParam} objects is returned.

For \code{checkList}, a character vector is returned containing all non-identical parameters across the list.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)
data$param

data <- windowCounts(bamFiles, width=100, filter=1, param=list(readParam(), readParam(minq=100)))
data$param
data[,1]$param
data[,2]$param

plist <- data$param
reformList(plist, restrict="chrA")
reformList(plist, minq=10)

checkList(plist)
checkList(reformList(plist, minq=10))
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}},
\code{\link{reform}}
}

\keyword{counting}
