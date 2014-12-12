# This is a testing stub that just provides a basic run-through of all the methods.
# The real tests are located in inst/tests along with a Bash script for execution
# and comparison. The separation into two folders is necessary to keep R CMD check
# running with reasonable speed.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))
both.files <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
pet.file <- system.file("exdata", "pet.bam", package="csaw")

# Checking data quality prior to counting.
head(correlateReads(both.files))
head(correlateReads(both.files, cross=FALSE))

getPETSizes(pet.file)

# Trying to count some single-end data.
data <- windowCounts(both.files, ext=100)
head(assay(data))
rowData(data)
data$total

data <- windowCounts(both.files, width=500, spacing=200)
head(assay(data))
rowData(data)
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
out <- windowCounts(pet.file, param=readParam(pet="both"), width=100, filter=1L)
assay(out)
out$totals
rowData(out)
out <- windowCounts(pet.file, param=readParam(pet="both", rescue.pairs=TRUE, rescue.ext=100), width=100, filter=1L)
assay(out)
out$totals
rowData(out)
out <- windowCounts(pet.file, param=readParam(pet="first"), width=100, filter=1L)
assay(out)
out$totals
rowData(out)

# Running some basic normalization.
data <- windowCounts(both.files, ext=100, param=readParam(minq=100, dedup=TRUE))
normalizeCounts(assay(data), lib.size=data$totals)
normalize(data)
normalizeCounts(assay(data), lib.size=data$totals, logratioTrim=.2)
normalizeCounts(assay(data), lib.size=data$totals, method="RLE")
head(normalizeCounts(assay(data), lib.size=data$totals, type="loess"))
head(normalize(data, type="loess"))
head(normalizeCounts(assay(data), lib.size=data$totals, type="loess", span=0.7))

# Assuming someone went around and pulled out some p-values for everybody.
set.seed(128145-19238)
nr <- nrow(data)
tabled <- data.frame(logFC=rnorm(nr), logCPM=rnorm(nr), PValue=rbeta(nr, 1, 2))
weighting <- rgamma(nr, 2, 1)

mergeWindows(rowData(data), -1)
mergeWindows(rowData(data), 100)
mergeWindows(rowData(data), 100, max.width=500)
merged <- mergeWindows(rowData(data), 100, sign=tabled$logFC > 0)
merged

head(combineTests(merged$id, tabled))
head(combineTests(merged$id, tabled, weight=weighting))

head(getBestTest(merged$id, tabled))
head(getBestTest(merged$id, tabled, weight=weighting))
head(getBestTest(merged$id, tabled, mode="logCPM"))

# Pulling out some diagnostics.
suppressPackageStartupMessages(require(org.Mm.eg.db))
suppressPackageStartupMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene))

current <- readRDS(system.file("exdata", "exrange.rds", package="csaw"))
output <- detailRanges(current, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db)
head(output$overlap)
head(output$right)
head(output$left)
     
output <- detailRanges(current, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, promoter=c(2000, 500))
head(output$overlap)
head(output$right)
head(output$left)

output <- detailRanges(current, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, dist=3000)
head(output$overlap)
head(output$right)
head(output$left)

# Pulling out some reads. 
extractReads(GRanges("chrA", IRanges(100, 500)), both.files[1])
extractReads(GRanges("chrA", IRanges(50, 100)), both.files[1])
extractReads(GRanges("chrA", IRanges(50, 100)), both.files[1], param=readParam(dedup=TRUE))
extractReads(GRanges("chrB", IRanges(50, 100)), pet.file, param=readParam(pet="both"))
extractReads(GRanges("chrB", IRanges(50, 100)), pet.file, param=readParam(pet="second"))

