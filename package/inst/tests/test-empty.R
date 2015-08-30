# Tests the behaviour of various functions when supplied with empty inputs.

suppressPackageStartupMessages(require(csaw))

empty <- data.frame(logFC=numeric(0), PValue=numeric(0), logCPM=numeric(0))
getBestTest(integer(0), empty)
getBestTest(integer(0), empty, by.pval=FALSE)
getBestOverlaps(Hits(), empty)

combineTests(integer(0), empty)
combineOverlaps(Hits(), empty)

upweightSummit(integer(0), integer(0))
summitOverlaps(Hits(), integer(0))

findMaxima(GRanges(), range=10, metric=numeric(0))

bamFile <- system.file("exdata", "rep1.bam", package="csaw")
profileSites(bamFile, GRanges(), range=20) # NA is correct, as average is undefined

checkBimodality(bamFile, GRanges())

clusterFDR(integer(0), 0.05) # NA is correct, as FDR is undefined.

suppressPackageStartupMessages(require(org.Mm.eg.db))
suppressPackageStartupMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene))
detailRanges(GRanges(), orgdb=org.Mm.eg.db, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene) 

out <- regionCounts(bamFile, GRanges())
out

getWidths(out)

asDGEList(out)

consolidateSizes(list(out), list(empty)) # No point testing behaviour on empty lists.
consolidateSizes(list(out), list(empty), region=GRanges())

reformList(list())
checkList(list())

makeExtVector(integer(0))

normOffsets(out) # 1 is correct, as calcNormFactors() just diverts to that.
normOffsets(out, type="loess")

scaledAverage(asDGEList(out))
scaledAverage(asDGEList(out), scale=numeric(0))

metadata(out)$spacing <- 50 # Converting to a window-based object.
chrs <- Rsamtools::scanBamHeader(bamFile)[[1]]$targets
seqinfo(rowRanges(out)) <- Seqinfo(names(chrs), chrs)
filterWindows(out, out, type="global")
filterWindows(out, out, type="local")
filterWindows(out, out, type="control", norm.fac=1) # TMM normalization fails for empty DGELists.
filterWindows(out, type="proportion")

