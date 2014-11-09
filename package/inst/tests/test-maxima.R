# This script tests the ability to identify local maxima from arbitary regions.

require(csaw)
source("simsam.R")

comp <- function(nregs, chromos, winsize, range) {
	reg.data <- generateWindows(chromos, nregs, winsize)
	reg.data$region <- reg.data$region[sample(nregs)]
	metric <- reg.data$table$logCPM

	# Getting max.
	obj <- SummarizedExperiment(matrix(0L, nregs, 1), rowData=reg.data$region,
		colData=DataFrame(row.names="A", whee=1))
	is.max <- findMaxima(obj, range=range, metric=metric)

	# Finding our own maxima.
	for (x in 1:length(nregs)) {
		new.reg <- reg.data$region[x]
		start(new.reg) <- start(new.reg) - range
		end(new.reg) <- end(new.reg) + range
		all.lap <- overlapsAny(reg.data$region, new.reg)
		
		check.max <- (metric[x] >= max(metric[all.lap]))
		if (!identical(check.max, is.max[x])) { stop("mismatch in max truths") }
	}

	return(sum(is.max))
}

#######################################################################################

set.seed(2394234)
comp(100, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=10, range=50)
comp(100, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=10, range=100)
comp(100, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=10, range=200)
comp(100, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=10, range=500)

comp(20, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=100, range=50)
comp(20, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=100, range=100)
comp(20, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=100, range=200)
comp(20, chromos=c(chrA=1000, chrB=2000, chrC=500), winsize=100, range=500)

comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=20, range=50)
comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=20, range=100)
comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=20, range=200)
comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=20, range=500)

comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=100, range=50)
comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=100, range=100)
comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=100, range=200)
comp(500, chromos=c(chrA=10000, chrB=20000, chrC=5000), winsize=100, range=500)

#######################################################################################

