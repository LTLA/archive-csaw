# This checks the overlap summarization functions, relative to the expected values.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))
source("simsam.R")

####################################################################################################

compcombine <- function(ranges, windows) {
	olap <- findOverlaps(ranges, windows)
	ns <- length(windows)
	tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))

	# Straight-up comparison to combineTests, after discarding all NA's.
	output <- combineOverlaps(olap, tab)
	refstats <- combineTests(queryHits(olap), tab[subjectHits(olap),])
	if (!identical(data.frame(output[!is.na(output$PValue),]), refstats)) { stop("mismatch in stats from near-identical calls!") }

	# Testing with weights.
	test.weight <- runif(ns)
	output <- combineOverlaps(olap, tab, i.weight=test.weight)
	refstats <- combineTests(queryHits(olap), tab[subjectHits(olap),], weight=test.weight[subjectHits(olap)])
	if (!identical(data.frame(output[!is.na(output$PValue),]), refstats)) { stop("mismatch in stats from near-identical calls!") }

	# More weight testing, where o.weight is constructed from the weight for each i.weight.
	output2 <- combineOverlaps(olap, tab, o.weight=test.weight[subjectHits(olap)])
	if (!identical(output, output2)) { stop("mismatch in stats from near-identical calls!") }

	return(head(output))
}

set.seed(34823)

chromos <- c(A=1000, B=2000)

regions <- generateWindows(chromos, 10, 500)
windows <- generateWindows(chromos, 100, 50)
compcombine(regions, windows)

regions <- generateWindows(chromos, 10, 500)
windows <- generateWindows(chromos, 200, 20)
compcombine(regions, windows)

regions <- generateWindows(chromos, 2, 500)
windows <- generateWindows(chromos, 100, 50)
compcombine(regions, windows)

####################################################################################################

compbest <- function(ranges, windows) {
	olap <- findOverlaps(ranges, windows)
	ns <- length(windows)
	tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))

	output <- getBestOverlaps(olap, tab)
	refstats <- getBestTest(queryHits(olap), tab[subjectHits(olap),])
	refstats$best <- subjectHits(olap)[refstats$best]
	if (!identical(data.frame(output[!is.na(output$PValue),]), refstats)) { stop("mismatch in stats from near-identical calls!") }

	# Testing with weights.
	test.weight <- runif(ns)
	output <- getBestOverlaps(olap, tab, i.weight=test.weight)
	refstats <- getBestTest(queryHits(olap), tab[subjectHits(olap),], weight=test.weight[subjectHits(olap)])
	refstats$best <- subjectHits(olap)[refstats$best]
	if (!identical(data.frame(output[!is.na(output$PValue),]), refstats)) { stop("mismatch in stats from near-identical calls!") }

	# More weight testing.
	output2 <- getBestOverlaps(olap, tab, o.weight=test.weight[subjectHits(olap)])
	if (!identical(output, output2)) { stop("mismatch in stats from near-identical calls!") }

	return(head(output))
}

set.seed(34823)

regions <- generateWindows(chromos, 10, 50)
windows <- generateWindows(chromos, 100, 50)
compbest(regions, windows)

regions <- generateWindows(chromos, 10, 50)
windows <- generateWindows(chromos, 200, 20)
compbest(regions, windows)

regions <- generateWindows(chromos, 2, 200)
windows <- generateWindows(chromos, 100, 50)
compbest(regions, windows)

####################################################################################################

compsummit <- function(ranges, windows) { 
	olap <- findOverlaps(ranges, windows)
	ns <- length(windows)
	tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
	output <- getBestOverlaps(olap, tab)

	# Checking summit calls.
	re.weight <- summitOverlaps(olap, output$best)
	best.win <- output$best[queryHits(olap)]
	is.summit <- !is.na(best.win) & best.win==subjectHits(olap)
	re.weight2a <- summitOverlaps(olap, o.summit=is.summit)
	re.weight2b <- summitOverlaps(olap, o.summit=which(is.summit))
	if (!identical(re.weight, re.weight2a) || !identical(re.weight, re.weight2b)) { stop("mismatch in weighting of summits") }

	isummits <- rbinom(ns, 1, 0.1)==1L
	re.weight3 <- summitOverlaps(olap, o.summit=isummits[subjectHits(olap)])
	re.weight4 <- summitOverlaps(olap, i.summit=isummits)
	if (!identical(re.weight3, re.weight4)) { stop("mismatch in weighting of summits") }

	# Checking the core upweightSummit machinery itself.
	by.region <- split(is.summit, queryHits(olap))	
	nu.weight <- sapply(by.region, FUN=function(x) {
		N <- length(x)
		output <- rep(1, N)
		output[x] <- N/sum(x)
		output	
	})
	if (!identical(re.weight, unlist(nu.weight, use.names=FALSE))) { stop("mismatch in manual weighting confirmation") }

	return(head(re.weight))
}

set.seed(34823)
	
regions <- generateWindows(chromos, 10, 50)
windows <- generateWindows(chromos, 20, 50)
compsummit(regions, windows)

regions <- generateWindows(chromos, 10, 50)
windows <- generateWindows(chromos, 20, 20)
compsummit(regions, windows)

regions <- generateWindows(chromos, 2, 200)
windows <- generateWindows(chromos, 10, 50)
compsummit(regions, windows)

####################################################################################################

