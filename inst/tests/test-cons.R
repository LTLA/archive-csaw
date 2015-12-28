# This tests the consolidateSizes function. As that function is basically written in fairly easy R, the
# test function below doesn't really do much except repeat the function itself.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))
source("simsam.R")
chromos <- c(chrA=12332, chrB=34892)

#########################################################################################################

compcons <- function(sizes, merge.args=list(tol=100), combine.args=list())  {
	data.list <- result.list <- list()
	for (s in 1:length(sizes)) {
		windows <- generateWindows(chromos, winsize=sizes[s], nwin=50)
		ns <- length(windows)
		result.list[[s]] <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
		
		counts <- matrix(0, ncol=1, nrow=ns)
		colnames(counts) <- "whee"
		data.list[[s]] <- SummarizedExperiment(assays=SimpleList(counts=counts), rowRanges=windows)
	}

	# Running consolidation.
	cons <- consolidateSizes(data.list, result.list, merge.args=merge.args, combine.args=combine.args)

	# Checking all regions.
	refbunch <- list()
	for (s in 1:length(sizes)) { 
		current <- rowRanges(data.list[[s]])
		comparator <- cons$region[cons$id[[s]]]
		if (any(seqnames(current)!=seqnames(comparator) | start(current) < start(comparator) | end(current) > end(comparator))) { 
				stop("merged regions do not include assigned windows")
		}
		refbunch[[s]] <- current
	}
	reference <- do.call(mergeWindows, c(regions=do.call(c, refbunch), merge.args))
	if (!identical(reference$region, cons$region)) { stop("merged regions don't match up with reference call") }
	refids <- unlist(cons$id)
	if (!identical(refids, reference$id)) { stop("merged ID's don't match up with reference call") }

	# Checking what happens with equi-weighting.
	new.weights <- list()
	for (s in 1:length(sizes)) { 
		counted <- tabulate(cons$id[[s]])
		new.weights[[s]] <- 1/counted[cons$id[[s]]]
	}	
	ref <- combineTests(reference$id, do.call(rbind, result.list), weight=unlist(new.weights))
	if (!identical(cons$table, ref)) { stop("mismatch in equiweighted results") }

	# Checking what happens without equi-weighting.		
	noe.cons <- consolidateSizes(data.list, result.list, merge.args=merge.args, combine.args=combine.args, equiweight=FALSE)
	ref <- combineTests(reference$id, do.call(rbind, result.list))
	if (!identical(noe.cons$table, ref)) { stop("mismatch in unweighted results") }
	
	return(cons$region)	
}

set.seed(2384)

compcons(c(100, 50))
compcons(c(100, 50), merge.args=list(tol=100, max.width=1000))
compcons(c(100, 50), merge.args=list(tol=1))
compcons(c(100, 50), combine.args=list(pval.col="PValue", fc.col="logFC"))

compcons(c(100, 500))
compcons(c(100, 500), merge.args=list(tol=100, max.width=1000))
compcons(c(100, 500), merge.args=list(tol=1))
try(compcons(c(100, 50), combine.args=list(weight=1)))

#########################################################################################################

compcons2 <- function(sizes, regions, overlap.args=list())  {
	data.list <- result.list <- list()
	for (s in 1:length(sizes)) {
		windows <- generateWindows(chromos, winsize=sizes[s], nwin=50)
		ns <- length(windows)
		result.list[[s]] <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
		counts <- matrix(0, ncol=1, nrow=ns)
		colnames(counts) <- "whee"
		data.list[[s]] <- SummarizedExperiment(assays=SimpleList(counts=counts), rowRanges=windows)
	}

	# Running consolidation.
	cons <- consolidateSizes(data.list, result.list, region=regions, overlap.args=overlap.args)

	# Checking all regions.
	refids <- list()
	for (s in 1:length(sizes)) { 
		current <- rowRanges(data.list[[s]])
		cur.lap <- do.call(findOverlaps, c(query=regions, subject=current, overlap.args))
		if (!identical(cons$id[[s]], cur.lap)) { stop("overlaps don't match up") } 
		refids[[s]] <- queryHits(cur.lap)
	}

	# Expanding out.
	neo.results <- list()
	for (s in 1:length(sizes)) { 
		neo.results[[s]] <- result.list[[s]][subjectHits(cons$id[[s]]),]
	}

	# Checking what happens with equi-weighting.
	new.weights <- list()
	for (s in 1:length(sizes)) { 
		counted <- tabulate(refids[[s]])
		new.weights[[s]] <- 1/counted[refids[[s]]]
	}	
	ref <- combineTests(unlist(refids), do.call(rbind, neo.results), weight=unlist(new.weights))
	if (!identical(cons$table[!is.na(cons$table$PValue),], ref)) { stop("mismatch in equiweighted results") }

	# Checking what happens without equi-weighting.		
	noe.cons <- consolidateSizes(data.list, result.list, region=regions, overlap.args=overlap.args, equiweight=FALSE)
	ref <- combineTests(unlist(refids), do.call(rbind, neo.results))
	if (!identical(noe.cons$table[!is.na(noe.cons$table$PValue),], ref)) { stop("mismatch in equiweighted results") }
	
	return(head(cons$table)	)
}

set.seed(2384)

regions <- generateWindows(chromos, winsize=1000, nwin=5)
compcons2(c(100, 50), regions)
compcons2(c(100, 50), regions, overlap.args=list(type="within"))
compcons2(c(100, 50), regions, overlap.args=list(minoverlap=10))

regions <- generateWindows(chromos, winsize=1000, nwin=5)
compcons2(c(10, 50), regions)
compcons2(c(10, 50), regions, overlap.args=list(type="within"))
compcons2(c(10, 50), regions, overlap.args=list(minoverlap=10))

regions <- generateWindows(chromos, winsize=1000, nwin=5)
compcons2(c(100, 500), regions)
compcons2(c(100, 500), regions, overlap.args=list(type="within"))
compcons2(c(100, 500), regions, overlap.args=list(minoverlap=10))

#########################################################################################################


