\name{wwhm}
\alias{wwhm}

\title{Window width at half maximum}
\description{Get the width of the window from the half-maximum of the coverage profile.}

\usage{
wwhm(profile, regions, ext=100, proportion=0.5, rlen=NULL)
}

\arguments{
\item{profile}{a numeric vector containing a coverage profile, as produced by \code{\link{profileSites}}}
\item{regions}{the \code{GRanges} object with which the profile was constructed}
\item{ext}{an integer scalar specifying the average fragment length for single-end data}
\item{proportion}{a numeric scalar specifying the proportion of the maximum coverage at which to compute the window width}
\item{rlen}{a numeric scalar or vector containing read lengths, if any \code{ext=NA}, i.e., fragments are unextended reads}
}

\details{
This function computes the ideal window size, based on the width of the peak in the coverage profile at the specified \code{proportion} of the maximum.
Obviously, the values of \code{regions} and \code{ext} should be the same as those used in \code{\link{profileSites}} (set \code{ext} to the median fragment length for paired-end data).
The \code{regions} should contain windows of a constant size.

Some subtleties are involved in obtaining the window width.
First, twice the average fragment length must be subtracted from the peak width, as the profile is constructed from (inferred) fragments.
The size of the viewpoints in \code{regions} must also be subtracted, to account for the inflated peak width when spatial resolution is lost after aggregation across many windows.

% Imagine if binding occured at 1 bp; random placement within the 50 bp bin would give a summit width of 50 bp.
}

\value{
An integer scalar is returned, specifying the ideal window width.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{profileSites}},
\code{\link{getWidths}}
}

\examples{
x <- dnorm(-200:200/100) # Mocking up a profile.
windows <- GRanges("chrA", IRanges(1, 50)) # Making up some windows.

wwhm(x, windows)
wwhm(x, windows, ext=50)
wwhm(x, windows, proportion=0.2)

# Need to set 'rlen' if ext=NA.
wwhm(x, windows, ext=NA, rlen=10)
}

\keyword{diagnostics}
