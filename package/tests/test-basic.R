# This is a testing stub that just provides a basic run-through of all the methods.
# The real tests are located in inst/tests along with a Bash script for execution
# and comparison. The separation into two folders is necessary to keep R CMD check
# running with reasonable speed.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))
both.files <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
pe.file <- system.file("exdata", "pet.bam", package="csaw")

# Checking data quality prior to counting.
head(correlateReads(both.files))
head(correlateReads(both.files, cross=FALSE))

getPESizes(pe.file)

# Trying to count some single-end data.
data <- windowCounts(both.files, ext=100)
head(assay(data))
rowRanges(data)
data$totals

data <- windowCounts(both.files, width=500, spacing=200)
head(assay(data))
rowRanges(data)
data$totals

data <- windowCounts(both.files, ext=100, param=readParam(minq=100))
data$totals
data <- windowCounts(both.files, ext=100, param=readParam(dedup=TRUE))
data$totals
data <- windowCounts(both.files, ext=100, param=readParam(discard=GRanges("chrA", IRanges(50, 500))))
data$totals
data <- windowCounts(both.files, ext=100, param=readParam(restrict="chrA"))
data$totals

# Trying to count some paired-end data.
out <- windowCounts(pe.file, param=readParam(pe="both"), width=100, filter=1L)
assay(out)
out$totals
rowRanges(out)
out <- windowCounts(pe.file, param=readParam(pe="both", rescue.ext=100), width=100, filter=1L)
assay(out)
out$totals
rowRanges(out)
out <- windowCounts(pe.file, param=readParam(pe="first"), width=100, filter=1L)
assay(out)
out$totals
rowRanges(out)

# Trying to convert it to a DGEList.
asDGEList(data)
asDGEList(data, lib.size=c(10, 10))
asDGEList(data, norm=c(1,2))

temp <- data 
temp$totals <- NULL 
asDGEList(temp) # Should spit out a warning.

# Running some basic normalization.
data <- windowCounts(both.files, ext=100, param=readParam(minq=100, dedup=TRUE))

normOffsets(assay(data), lib.size=data$totals)
normOffsets(assay(data), lib.size=data$totals, logratioTrim=.2)
normOffsets(assay(data), lib.size=data$totals, method="RLE")
normOffsets(data)
normOffsets(data, logratioTrim=0.1)
normOffsets(data, method="upperquartile")

head(normOffsets(assay(data), lib.size=data$totals, type="loess"))
head(normOffsets(assay(data), lib.size=data$totals, type="loess", span=0.7))
head(normOffsets(data, type="loess"))
head(normOffsets(data, type="loess", span=0.5))

# Assuming someone went around and pulled out some p-values for everybody.
set.seed(128145-19238)
nr <- nrow(data)
tabled <- data.frame(logFC=rnorm(nr), logCPM=rnorm(nr), PValue=rbeta(nr, 1, 2))
weighting <- rgamma(nr, 2, 1)

mergeWindows(rowRanges(data), -1)
mergeWindows(rowRanges(data), 100)
mergeWindows(rowRanges(data), 100, max.width=500)
merged <- mergeWindows(rowRanges(data), 100, sign=tabled$logFC > 0)
merged

head(combineTests(merged$id, tabled))
head(combineTests(merged$id, tabled, weight=weighting))

head(getBestTest(merged$id, tabled))
head(getBestTest(merged$id, tabled, weight=weighting))
head(getBestTest(merged$id, tabled, by.pval=FALSE))

# Pulling out some diagnostics.
suppressPackageStartupMessages(require(org.Mm.eg.db))
suppressPackageStartupMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene))

current <- readRDS(system.file("exdata", "exrange.rds", package="csaw"))
output <- detailRanges(current, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db)
head(output$overlap)
head(output$right)
head(output$left)

# Pulling out some reads. 
extractReads(both.files[1], GRanges("chrA", IRanges(100, 500)))
extractReads(both.files[1], GRanges("chrA", IRanges(50, 100)))
extractReads(both.files[1], GRanges("chrA", IRanges(50, 100)), param=readParam(dedup=TRUE))
extractReads(pe.file, GRanges("chrB", IRanges(50, 100)), param=readParam(pe="both"))
extractReads(pe.file, GRanges("chrB", IRanges(50, 100)), param=readParam(pe="second"))

