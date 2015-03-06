\name{csawnews}
\title{csaw News}
\encoding{UTF-8}

\section{Version 1.1.31}{\itemize{
\item 
Added getWidths(), scaledAverage() and filterWindows(), to facilitate comparison of abundances during filtering.

\item 
Added findMaxima() to identify locally maximal windows.

\item 
Added profileSites() to examine the coverage profile around specified regions, with wwhm() to guess the ideal window size from the profile.

\item
Changed default window width in windowCounts() to 50 bp, default filter to 10.

\item 
Changed all instances of \code{pet} to \code{pe} in read parameter specification, and renamed getPETSizes() to getPESizes().

\item 
Removed the redundant rescue.pairs parameter in readParam().

\item
Added fast.pe parameter in readParam(), for fast paired-end data extraction. Added dumpPE() to pre-process paired-end BAM files for fast downstream extraction.

\item 
Added support for custom column specification in getBestTest(), combineTests().

\item 
Added consolidateSizes(), to consolidate DB results from multiple window sizes.

\item
Added upweightSummits(), to compute weights favouring high-abundance windows during Simes.

\item
Added support for custom key/name specification in detailRanges() for non-human/mouse systems.

\item 
Added support for strand-specific read extraction in readParam(), strand-specific counting via strandedCounts().

\item 
Added strand-awareness to mergeWindows(). Added protection against stranded input regions in extractReads(), detailRanges().

\item
Changed algorithm for splitting of large peaks in mergeWindows().

\item
Stored counting parameters in exptData for windowCounts(), regionCounts().

\item
Fixed small inaccuracies with continuity correction addition in normalizeCounts() for NB-loess.

\item 
Switched to fragment midpoint for binning of paired-end data in windowCounts().

\item 
Added support for lists of library-specific readParam objects in windowCounts(), regionCounts(), correlateReads().

\item 
Added makeExtVector(), to support variable read extension lengths between libraries in windowCounts(), regionCounts().

\item
Added support for read extension within extractReads().

\item 
Updated the user's guide to reflect new and modified functions.

\item 
Cleaned up code in inst/tests for modified functions, added new tests for new functions.
}}

\section{Version 1.0.0}{\itemize{
\item 
New package \code{csaw}, for de novo detection of differential binding in ChIP-seq data.
}}
