\name{Parameter list methods}
\alias{paramList}
\alias{reformList}
\alias{checkList}

\title{Get or modify readParam lists}
\description{Extract and modify lists of library-specific readParam objects}

\usage{
paramList(data)
reformList(paramlist, ...)
checkList(paramlist)
}

\arguments{
  \item{data}{a \code{SummarizedExperiment} object, like that produced by \code{\link{windowCounts}}}
  \item{paramlist}{a list of \code{readParam} objects}
  \item{...}{other arguments, to be passed to \code{\link{reform}}}
}

\details{
The \code{paramList} function extracts a list of \code{readParam} objects from a \code{SummarizedExperiment}.
Each object represents the parameter settings for each library used in counting by \code{\link{windowCounts}} or \code{\link{regionCounts}}.

The \code{reformList} function sets the parameters specified in \code{...} across all elements of \code{paramlist}.
This standardizes all of the internal \code{readParam} objects, e.g., for the \code{restrict} parameter.
If \code{paramlist} is a single \code{readParam} object, then \code{reformList} is equivalent to \code{\link{reform}}.

The \code{checkList} function checks which parameter settings are identical between elements of \code{paramlist}.
This is useful when some parameters must be variable, while others must be the same.
}

\value{
For \code{paramList} and \code{reformList}, a list of library-specific \code{readParam} objects is returned.

For \code{checkList}, a character vector is returned containing all non-identical parameters across the list.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)
paramList(data)

data <- windowCounts(bamFiles, width=100, filter=1, param=list(readParam(), readParam(minq=100)))
paramList(data)
paramList(data[,1])
paramList(data[,2])

plist <- paramList(data)
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

