\name{csawUsersGuide}
\alias{csawUsersGuide}

\title{View csaw user's guide}
\description{Finds the location of the user's guide and opens it for viewing.}

\usage{
csawUsersGuide(view=TRUE)
}

\arguments{
  \item{view}{logical scalar specifying whether the document should be opened}
}

\value{A character string giving the file location. If \code{view=TRUE}, the
system's default PDF document reader is started and the user's guide is 
opened.}

\details{
Calling \code{vignette("csaw")} will yield a short vignette that contains
little information. Instead, the user guide can be obtained with this function.
The user's guide is not a true vignette as it is not generated using
\code{\link{Sweave}} when the package is built. This is due to the
time-consuming nature of the code when run on realistic case studies.

For non-Windows operating systems, the PDF viewer is taken from
\code{Sys.getenv("R_PDFVIEWER")}. This can be changed to \code{x} by using
\code{Sys.putenv(R_PDFVIEWER=x)}. For Windows, the default viewer will be
selected to open the file. Note that for Windows, the user's guide can
also be accessed from the `Vignettes' drop-down menu in the R GUI.
}

\seealso{
\code{\link{system}}
}

\author{Aaron Lun}

\examples{
# To get the location:
csawUsersGuide(view=FALSE)
# To open in pdf viewer:
\dontrun{csawUsersGuide()}
}

\keyword{documentation}
