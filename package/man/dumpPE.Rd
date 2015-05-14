\name{dumpPE}
\alias{dumpPE}

\title{Dump paired-end data to file}
\description{Extract proper pairs from a BAM file, and dump fragment intervals into another BAM file.}

\usage{
dumpPE(bam.file, prefix, param=readParam(pe="both"), overwrite=FALSE)
}

\arguments{
\item{bam.file}{a character string containing the path to a paired-end BAM file}
\item{prefix}{a character string containing the prefix to an output BAM file}
\item{param}{a \code{readParam} object}
\item{overwrite}{a logical scalar indicating whether to overwrite the existing file starting with \code{prefix}}
}

\value{
A sorted and indexed BAM file is produced at the specified location.
A character string containing the full name of the output BAM file is silently returned.
}

\details{
This function extracts proper pairs from \code{bam.file} according to the settings in \code{param}.
It then generates another output BAM file that stores fragment information for each proper pair.
Each alignment entry represents the forward-stranded read of a valid fragment.
Fragment data can be extracted as the read position and insert size.

The idea is to generate a pre-filtered BAM file by specifying the appropriate settings in \code{param}.
The output BAM file can then be efficiently analyzed with, e.g., \code{\link{windowCounts}} with \code{fast.pe=TRUE}.
This avoids the need to load and match read names in every function.
Note that all alignment-specific information is lost, e.g., MAPQ scores, duplicate information, CIGAR strings.

Only proper pairs with sizes less than \code{param$max.frag} will be retained in the output file.
If \code{max.frag} is increased in downstream calls on the output file, the results may not be correct as proper pairs with sizes between the two settings will not be considered.
Users are recommended to keep \code{max.frag} constant throughout the pipeline.
}

\seealso{
\code{\link{readParam}},
\code{\link{windowCounts}}
}

\examples{
# Loading PE data.
bamFile <- system.file("exdata", "pet.bam", package="csaw")

xparam <- readParam(pe="both")
out <- windowCounts(bamFile, param=xparam, filter=1)

outBam <- dumpPE(bamFile, "whee", param=xparam)
out2 <- windowCounts(outBam, param=reform(xparam, fast.pe=TRUE), filter=1)

stopifnot(identical(assay(out), assay(out2)))
stopifnot(identical(out$totals, out2$totals))

# Comparing it to a more complex scenario.
xparam <- readParam(pe="both", rescue.ext=200)
out <- windowCounts(bamFile, param=xparam, filter=1)

outBam <- dumpPE(bamFile, "whee", param=xparam, overwrite=TRUE)
out2 <- windowCounts(outBam, param=reform(xparam, fast.pe=TRUE), filter=1)

stopifnot(identical(assay(out), assay(out2)))
stopifnot(identical(out$totals, out2$totals))

# Looking at extracted reads.
last.reg <- rowRanges(out)[6]
reg1 <- extractReads(bamFile, last.reg, param=xparam)
reg2 <- extractReads(outBam, last.reg, param=reform(xparam, fast.pe=TRUE))

stopifnot(identical(sort(reg1), sort(reg2)))

\dontshow{
unlink("whee.bam")
unlink("whee.bam.bai")
}
}

\keyword{counting}
