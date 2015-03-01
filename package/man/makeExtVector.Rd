\name{makeExtVector}
\alias{makeExtVector}
\alias{[.extVector}

\title{Make an extVector object}
\description{Construct an extVector to hold read extension parameters.}

\usage{
makeExtVector(ext, final.ext=NULL)
}

\arguments{
\item{ext}{an integer scalar or vector containing average fragment lengths for read extension in one or more libraries}
\item{final.ext}{an integer scalar, specifying the length to which all fragments are scaled}
}

\value{
An \code{extVector} object containing the values of \code{ext}, with the value of \code{final.ext} stored in the attributes.
}

\details{
This function returns an \code{extVector} object, which is just an augmented integer vector with the \code{final.ext} attribute.
The aim of using a specialized class is to preserve attributes during subsetting.
Thus, functions like \code{\link{windowCounts}} or \code{\link{regionCounts}} can be run with different combinations of libraries, and \code{final.ext} will be maintained in each call.

The values in \code{ext} will ultimately be used for directional extension, to infer the fragment corresponding to each read.
Different fragment lengths will be used in different libraries if a vector is supplied for \code{ext}.
If \code{final.ext} is \code{NULL}, it is defined by \code{makeExtVector} as the mean of all \code{ext}.

Lengths of all fragments will be scaled to the same \code{final.ext} value in each library.
If \code{final.ext} is \code{NA}, no scaling is performed.
See \code{\link{windowCounts}} for more details. 
}

\seealso{
\code{\link{windowCounts}},
\code{\link{regionCounts}}
}

\examples{
out <- makeExtVector(1:5*10, NA)
out
out[1:2]

out <- makeExtVector(1:5*10, NULL)
out
out[1:2] # Doesn't recalculate, which is correct.
}

\author{Aaron Lun}

\keyword{counting}
