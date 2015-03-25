\name{upweightSummit}
\alias{upweightSummit}

\title{Upweight summits}
\description{Upweight the highest-abudance window(s) in a cluster.}

\usage{
upweightSummit(ids, summits)
}

\arguments{
\item{ids}{an integer vector of cluster IDs}
\item{summits}{a logical vector indicating whether each window is a summit, or an integer vector containing the indices of summit windows}
}

\value{
A numeric vector of weights, where the highest-abundance window in each cluster is assigned a greater weight.
}

\details{
This function computes weights for each window in a cluster, where the highest-abundance windows are upweighted.
These weights are intended for use in \code{\link{combineTests}}, such that the summits of a cluster have a greater influence on the combined p-value.
This is more graduated than simply using the summits alone, as potential DB between summits can still be detected.
Summits can be obtained through \code{\link{findMaxima}} or by running \code{\link{getBestTest}} with \code{by.pval=FALSE}.

The exact value of the weight is arbitrary.
Greater weight represents a stronger belief that DB occurs at the most abundant window.
Here, the weighting scheme is designed such that the maximum Simes correction is not more than twice that without weighting.
It will also be no more than twice that from applying Simes' method on the summits alone.
This (restrained) conservativeness is an acceptable cost for considering DB events elsewhere in the cluster, while still focusing on the most abundant site.

% To wit; let 'w' be the weight of each summit (equal weights), S be the number of summits and N be the number of windows in a cluster.
% If one of the summits is the most DB, then the extra penalty is (N+S*w-S)/(S*w), i.e., the multiplicity correction compared to just using the p-values of the summits.
% If another window is most DB, then the extra penalty is (N+S*w-S)/N, i.e., compared to unweighted Simes.
% These two penalties can be made equal by setting w=N/S, which means the penalty is below 2 (how much depends on the amount of S relative to N).
}

\seealso{
\code{\link{combineTests}},
\code{\link{findMaxima}},
\code{\link{getBestTest}}
}

\references{
Benjamini Y and Hochberg Y (1997). Multiple hypotheses testing with weights. \emph{Scand. J. Stat.} 24, 407-418.
}

\author{
Aaron Lun
}

\examples{
nwin <- 20
set.seed(20)
ids <- sample(5, nwin, replace=TRUE)	
summits <- sample(5, nwin, replace=TRUE)==1L
weights <- upweightSummit(ids, summits)

# Checking that the summit is upweighted in each cluster.
split(data.frame(summits, weights), ids)
}

\keyword{testing}
