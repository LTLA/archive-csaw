\name{csawnews}
\title{csaw News}
\encoding{UTF-8}

\section{Version 1.1.11}{\itemize{
\item 
Added getWidths() and scaledAverage(), to facilitate comparison of abundances during filtering.

\item 
Added findMaxima() to identify locally maximal windows.

\item 
Added profileSites() to examine the coverage profile around specified regions.

\item
Changed default width in windowCounts() to 50.

\item
Changed algorithm for splitting of large peaks in mergeWindows().

\item
Stored counting parameters in exptData for windowCounts(), regionCounts().

\item
Fixed small inaccuracies with continuity correction addition in normalizeCounts() for NB-loess.

\item 
Switched to fragment midpoint for binning of paired-end data in windowCounts().

\item 
Updated the user's guide to reflect new and modified functions.

\item 
Cleaned up code in inst/tests for modified functions, added new tests for new functions.
}}

\section{Version 1.0.0}{\itemize{
\item 
New package \code{csaw}, for de novo detection of differential binding in ChIP-seq data.
}}
