\name{filterWindows}
\alias{filterWindows}

\title{Filtering methods for SummarizedExperiment objects}
\description{Convenience function to compute filter statistics for windows, based on proportions or using enrichment over background.}

\usage{
filterWindows(data, background, type="global", prior.count=2, len=NULL) 
}

\arguments{
\item{data}{a \code{SummarizedExperiment} object containing window- or bin-level counts}
\item{background}{another \code{SummarizedExperiment} object, containing counts for background regions when \code{type!="proportion"}}
\item{type}{a character string specifying the type of filtering to perform; can be any of \code{c("global", "local", "control", "proportion")}}
\item{prior.count}{a numeric scalar, specifying the prior count to use in \code{\link{aveLogCPM}}}
\item{len}{a numeric scalar containing the median fragment length for paired-end data}
}

\value{
A list is returned with \code{abundances}, the average abundance of each entry in \code{data}; \code{filter}, the filter statistic for the given \code{type}; and, for \code{type!="proportion"}, \code{back.abundances}, the average abundance of each entry in \code{background}.
}

\details{
Proportion-based filtering supposes that a certain percentage of the genome is genuinely bound.
If \code{type="proportion"}, the filter statistic is defined as the ratio of the rank to the total number of windows.
Rank is in ascending order, i.e., higher abundance windows have higher ratios.
Windows are retained that have rank ratios above a threshold, e.g., 0.99 if 1\% of the genome is assumed to be bound.

All other values of \code{type} will perform background-based filtering, where abundances of the windows are compared to those of putative background regions.
The filter statistic are generally defined as the difference between window and background abundances, i.e., the log-fold increase in the counts.
Windows can be filtered to retain those with large filter statistics, to select those that are more likely to contain genuine binding sites.
The differences between the methods center around how the background abundances are obtained for each window.

If \code{type="global"}, the median average abundance across the genome is used as a global estimate of the background abundance.
This should be used when \code{background} contains unfiltered counts for large (2 - 10 kbp) genomic bins, from which the background abundance can be computed.
The filter statistic for each window is defined as the difference between the window abundance and the global background.
If \code{background} is not supplied, the background abundance is directly computed from entries in \code{data}.

% Note that empty windows are automatically removed by windowCounts.
% However, if you're using large bins, these regions are likely to be systematically empty, e.g. large repeats like telomeres, centromeres.
% I'm happy to consider them as being irrelevant parts of the genome that we can ignore, as they don't help in estimating the background (as they're empty).

If \code{type="local"}, the counts of each row in \code{data} are subtracted from those of the corresponding row in \code{background}.
The average abundance of the remaining counts is computed and used as the background abundance.
The filter statistic is defined by subtracting the background abundance from the corresponding window abundance for each row.
This is designed to be used when \code{background} contains counts for expanded windows, to determine the local background estimate.

If \code{type="control"}, the background abundance is defined as the average abundance of each row in \code{background}.
The filter statistic is defined as the difference between the average abundance of each row in \code{data} and that of the corresponding row in \code{background}.
This is designed to be used when \code{background} contains read counts for each window in the control sample(s).
Unlike \code{type="local"}, there is no subtraction of the counts in \code{background} prior to computing the average abundance.
}

\section{Additional details}{
Proportion and global background filtering are dependent on the total number of windows/bins in the genome.
However, empty windows or bins are automatically discarded in \code{\link{windowCounts}} (exacerbated if \code{filter} is set above unity).
This will result in underestimation of the rank or overestimation of the global background.
To avoid this, the total number of windows or bins is inferred from the spacing.

For background-based methods, the abundances of large bins or regions in \code{background} must be rescaled for comparison to those of smaller windows 
- see \code{\link{getWidths}} and \code{\link{scaledAverage}} for more details.
In particular, the effective width of the window is often larger than \code{width}, due to the counting of fragments rather than reads.
The fragment length is generally extracted from \code{data}.
If not (e.g., for paired-end analyses), specification of \code{len} may be necessary.

The \code{prior.count} protects against inflated log-fold increases when the background counts are near zero.
A low prior is sufficient if \code{background} has large counts, which is usually the case for wide regions.
Otherwise, \code{prior.count} should be increased to a larger value like 5.
This may be necessary in \code{type="control"}, where \code{background} contains counts for small windows in the control sample.
}

\seealso{
\code{\link{windowCounts}},
\code{\link{aveLogCPM}},
\code{\link{getWidths}},
\code{\link{scaledAverage}}
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, filter=1)

# Proportion-based (keeping top 1%)
stats <- filterWindows(data, type="proportion")
head(stats$filter)
keep <- stats$filter > 0.99 
new.data <- data[keep,]

# Global background-based (keeping fold-change above 3).
background <- windowCounts(bamFiles, bin=TRUE, width=300)
stats <- filterWindows(data, background, type="global")
head(stats$filter)
keep <- stats$filter > log2(3)

# Local background-based.
locality <- regionCounts(bamFiles, resize(rowData(data), fix="center", 300))
stats <- filterWindows(data, locality, type="local")
head(stats$filter)
keep <- stats$filter > log2(3)

# Control-based (pretend "rep.2" is a control library).
stats <- filterWindows(data[,1], data[,2], type="control", prior.count=5)
head(stats$filter)
keep <- stats$filter > log2(3)
}

\keyword{filtering}
