\name{filterWindows}
\alias{filterWindows}

\title{Filtering methods for SummarizedExperiment objects}
\description{Convenience function to compute filter statistics for windows, based on proportions or using enrichment over background.}

\usage{
filterWindows(data, background, type="global", prior.count=2) 
}

\arguments{
\item{data}{a \code{SummarizedExperiment} object containing window- or bin-level counts}
\item{background}{another \code{SummarizedExperiment} object, containing counts for background regions when \code{type!="proportion"}}
\item{type}{a character string specifying the type of filtering to perform; can be any of \code{c("global", "local", "control", "proportion")}}
\item{prior.count}{a numeric scalar, specifying the prior count to use in \code{\link{aveLogCPM}}}
}

\value{
A list is returned with \code{abundances}, the average abundance of each entry in \code{data}; \code{filter}, the filter statistic for the given \code{type}; and, for \code{type!="proportion"}, \code{back.abundances}, the average abundance of each entry in \code{background}.
}

\details{
Proportion-based filtering supposes that a certain percentage of the genome is genuinely bound.
If \code{type="proportion"}, the filter statistic is defined as the ratio of the rank to the total number of windows.
Rank is in ascending order, i.e., higher abundance windows have higher ratios.
Determining the total number requires some care, as the majority of the genome is already discarded by pre-filtering in \code{\link{windowCounts}}.
Windows are retained that have rank ratios above a threshold, e.g., 0.99 if 1\% of the genome is assumed to be bound.

Background-based filtering compares the window abundances to those of putative background regions.
If \code{type="global"}, the median average abundance of all entries in \code{background} is used as a global estimate of the background abundance.
This is designed to be used when \code{background} contains counts for large (2000 - 10000 kbp) genomic bins.
The filter statistic is defined as the difference between the window abundance and the global background, i.e., the log-fold increase over the background.

If \code{type="local"}, the counts of each row in \code{data} are subtracted from those of the corresponding row in \code{background}.
The average abundance of the remaining counts is computed and used as the background abundance.
The filter statistic is defined by subtracting the background abundance from the corresponding window abundance for each row.
This is designed to be used when \code{background} contains counts for expanded windows, to determine the loal background estimate.

If \code{type="control"}, the filter statistic is defined by subtracting the average abundance of each row in \code{background} from the average abundance of the corresponding row in \code{data}.
This is designed to be used when \code{background} contains read counts for each window in the control sample(s).

For all background-based methods, windows can be filtered to retain those with large filter statistics.
This selects windows with strong enrichment over background, that are more likely to contain genuine binding sites.
In all cases, the abundances of large bins or regions in \code{background} must be rescaled for comparison to those of smaller windows.
}

\seealso{
\code{\link{windowCounts}},
\code{\link{aveLogCPM}}
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, filter=1)

# Proportion-based (keeping top 1%)
stats <- filterWindows(data, type="proportion")
head(stats$filter)
keep <- stats$filter > 0.99 
new.data <- data[keep,]

# Global background-based (keeping fold-change above 5).
background <- windowCounts(bamFiles, bin=TRUE, width=300)
stats <- filterWindows(data, background, type="global")
head(stats$filter)
keep <- stats$filter > log2(5)

# Local background-based.
locality <- regionCounts(bamFiles, resize(rowData(data), fix="center", 300))
stats <- filterWindows(data, locality, type="local")
head(stats$filter)
keep <- stats$filter > log2(5)

# Control-based (pretend "rep.2" is a control library).
stats <- filterWindows(data[,1], data[,2], type="control")
head(stats$filter)
keep <- stats$filter > log2(5)
}

