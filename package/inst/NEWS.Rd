\name{csawnews}
\title{csaw News}
\encoding{UTF-8}

\section{Version 1.3.1}{\itemize{
\item
Added clusterFDR() function to compute the FDR for clusters of DB windows.

\item
Modified normalize(), asDGEList() to allow manual specification of library sizes.

\item 
Modified default parameter specification in strandedCounts(), to avoid errors.

\item
Slight updates to documentation, user's guide.
}}

\section{Version 1.2.0}{\itemize{
\item 
Added getWidths(), scaledAverage() and filterWindows(), to facilitate comparison of abundances during filtering.

\item 
Added findMaxima() to identify locally maximal windows from range data.

\item 
Added profileSites() to examine the coverage profile around specified regions, with wwhm() to guess the ideal window size from the profile.

\item
Changed default window width in windowCounts() to 50 bp, default filter to a fixed count of 10. Also, filter=0 is honored when bin=TRUE.

\item
Switched from the depracated rowData to rowRanges for all manipulations of SummarizedExperiment.

\item 
Changed all instances of `pet' to `pe' in read parameter specification, and renamed getPETSizes() to getPESizes().

\item 
Removed the redundant rescue.pairs parameter in readParam().

\item
Added fast.pe parameter in readParam(), for fast paired-end data extraction. Added dumpPE() to pre-process paired-end BAM files for fast downstream extraction.

\item 
Added support for custom column specification in getBestTest(), combineTests().

\item
Switched from reporting average log-FC to numbers of up/down windows in combineTests().

\item
Allowed getBestTest() to return all fields associated with the best window in the output table.

\item
Added upweightSummits() to compute weights favouring high-abundance windows.

\item 
Added combineOverlaps(), getBestOverlaps() and summitOverlaps() wrapper functions for processing of Hits objects.

\item 
Added consolidateSizes(), to consolidate DB results from multiple window sizes.

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
Added sra2bam.sh in inst/doc to reproducibly generate BAM files to run UG examples.

\item 
Cleaned up code in inst/tests for modified functions, added new tests for new functions.
}}

\section{Version 1.0.0}{\itemize{
\item 
New package \code{csaw}, for de novo detection of differential binding in ChIP-seq data.
}}
