# This is a testing stub that just provides a basic run-through of all the methods.
# The real tests are located in inst/tests along with a Bash script for execution
# and comparison. The separation into two folders is necessary to keep R CMD check
# running with reasonable speed.

suppressPackageStartupMessages(require(csaw))
both.files <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
pet.file <- system.file("exdata", "pet.bam", package="csaw")

# Checking data quality prior to counting.
head(correlateReads(both.files))
head(correlateReads(both.files, cross=FALSE))

getPETSizes(pet.file)

# Trying to count some single-end data.
data <- windowCounts(both.files, ext=100)
head(data$counts)
data$region
data$total

data <- windowCounts(both.files, width=500, spacing=200)
head(data$counts)
data$region
data$totals

data <- windowCounts(both.files, ext=100, minq=100)
data$totals
data <- windowCounts(both.files, ext=100, dedup=TRUE)
data$totals
data <- windowCounts(both.files, ext=100, discard=GRanges("chrA", IRanges(50, 500)))
data$totals
data <- windowCounts(both.files, ext=100, restrict="chrA")
data$totals

param <- list(bam.files=both.files, ext=100, minq=100, dedup=TRUE)
data <- windowCounts(both.files, ext=100, minq=100, dedup=TRUE)
identical(countWindows(param), data) 
identical(countWindows(param, ext=200), windowCounts(both.files, ext=200, minq=100, dedup=TRUE))
identical(countWindows(param, ext=50, minq=0, dedup=FALSE), windowCounts(both.files, ext=50, minq=0, dedup=FALSE))

# Trying to count some paired-end data.
windowCounts(pet.file, pet="both", width=100, filter=1L)
windowCounts(pet.file, pet="both", width=100, filter=1L, rescue.pairs=TRUE)
windowCounts(pet.file, pet="first", width=100, filter=1L)

# Running some basic normalization.
normalizeChIP(data$counts, lib.size=data$totals)
normalizeChIP(data$counts, lib.size=data$totals, logratioTrim=.2)
normalizeChIP(data$counts, lib.size=data$totals, method="RLE")
head(normalizeChIP(data$counts, lib.size=data$totals, type="loess"))
head(normalizeChIP(data$counts, lib.size=data$totals, type="loess", span=0.7))

# Assuming someone went around and pulled out some p-values for everybody.
set.seed(128145-19238)
nr <- nrow(data$counts)
tabled <- data.frame(logFC=rnorm(nr), logCPM=rnorm(nr), PValue=rbeta(nr, 1, 2))
weighting <- rgamma(nr, 2, 1)

mergeWindows(data$region, 10)
mergeWindows(data$region, 100)
mergeWindows(data$region, 100, max.width=500)
merged <- mergeWindows(data$region, 100, sign=tabled$logFC > 0)
merged

head(combineTests(merged$id, tabled))
head(combineTests(merged$id, tabled, weight=weighting))

head(getBestTest(merged$id, tabled))
head(getBestTest(merged$id, tabled, weight=weighting))
head(getBestTest(merged$id, tabled, mode="logCPM"))

# Pulling out some diagnostics.
require(org.Mm.eg.db)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)

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

# Here's a pretty plot.
pdf("tracktest.pdf")
plotRegion(GRanges("chrA", IRanges(100, 500)), both.files[1])
plotRegion(GRanges("chrA", IRanges(50, 100)), both.files[1])
plotRegion(GRanges("chrB", IRanges(50, 100)), pet.file)
dev.off()

unlink("tracktest.pdf")
