\name{maximizeCcf}
\alias{maximizeCcf}

\title{Find the delay at the maximal CCF}
\description{Estimate the average fragment length by maximizing the cross-correlations.}

\usage{
maximizeCcf(profile, ignore=100)
}

\arguments{
\item{profile}{a numeric vector containing a coverage profile, as produced by \code{\link{correlateReads}}}
\item{ignore}{an integer scalar specifying the distances to ignore}
}

\details{
This function identifies the delay distance at which the cross-correlations are maximized.
This distance can then be used as an estimate of the average fragment length, for use in directional extension during read counting. 

In some datasets, identification of the maxima is confounded by a phantom peak at the read length.
This can be overcome by ignoring the first \code{ignore} delay distances, such that the distance corresponding to the true peak is used.

Obviously, this only works in TF experiments with moderate to strong enrichment, where a strong peak in the CCF profile is present. 
The function may not perform sensibly in the presence of noisy profiles containing multiple local maxima.
}

\value{
The average fragment length is returned as an integer scalar.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{correlateReads}}
}

\references{
Landt SG, Marinov GK, Kundaje A, et al. (2012).
ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia.
\emph{Genome Res.} 22, 1813-31.
}

\examples{
x <- dnorm(-200:200/100) # Mocking up a profile.
maximizeCcf(x)

x2 <- x + dnorm(-50:250/10) # Adding a phantom peak
maximizeCcf(x2)
maximizeCcf(x2, ignore=0)
}

\keyword{diagnostics}
