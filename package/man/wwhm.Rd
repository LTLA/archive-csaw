\name{wwhm}
\alias{wwhm}

\title{Window width at half maximum}
\description{Get the width of the window from the half-maximum of the coverage profile.}

\usage{
wwhm(profile, regions, ext=100, param=readParam(), proportion=0.5, len=NULL)
}

\arguments{
\item{profile}{a numeric vector containing a coverage profile, as produced by \code{\link{profileSites}}}
\item{regions}{the \code{GRanges} object with which the profile was constructed}
\item{ext}{an integer scalar specifying the average fragment length for single-end data}
\item{param}{a \code{readParam} object containing read extraction parameters, or a list of such objects}
\item{proportion}{a numeric scalar specifying the proportion of the maximum coverage at which to compute the window width}
\item{len}{a numeric scalar or vector containing average fragment lengths, to be used in \code{\link{getWidths}} in special cases}
}

\details{
This function computes the ideal window size, based on the width of the peak in the coverage profile at the specified \code{proportion} of the maximum.
Obviously, the values of \code{regions}, \code{ext} and \code{param} should be the same as those used in \code{\link{profileSites}}.
The \code{regions} should contain windows of a constant size.

Some subtleties are involved in obtaining the window width.
First, twice the average fragment length must be subtracted from the peak width, as the profile is constructed from (inferred) fragments.
The size of the viewpoints in \code{regions} must also be subtracted, to account for the inflated peak width when spatial resolution is lost.

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
}

\keyword{diagnostics}
