# This tests the post-hoc clustering methods.

suppressPackageStartupMessages(require(csaw))
compCFDR <- function(ids, threshold, weights) {
    if (is.null(weights)) { 
        obs.sizes <- table(ids) 
    } else {
        obs.sizes <- sapply(split(weights, ids), FUN=sum)
    }
    obs.sizes <- sort(obs.sizes)
    num.fp <- sum(cumsum(obs.sizes) <= sum(obs.sizes) * threshold)
    ref.fdr <- num.fp/length(obs.sizes)
    test.fdr <- clusterFDR(ids, threshold, weight=weights)
    stopifnot(all.equal(test.fdr, ref.fdr))
    return(test.fdr)
}

set.seed(100)
ids <- rnbinom(100, mu=5, size=10)
compCFDR(ids, 0.05, NULL)
compCFDR(ids, 0.1, NULL)
compCFDR(ids, 0.05, runif(100))

ids <- rnbinom(100, mu=10, size=1)
compCFDR(ids, 0.05, NULL)
compCFDR(ids, 0.1, NULL)
compCFDR(ids, 0.05, runif(100))

ids <- rnbinom(100, mu=100, size=20)
compCFDR(ids, 0.05, NULL)
compCFDR(ids, 0.1, NULL)
compCFDR(ids, 0.05, runif(100))

###################################################
# Checking that the weighted p-value methods are correct.

set.seed(1000)
pvals <- rbeta(1000, 1, 20)
stopifnot(all.equal(p.adjust(pvals, method="BH"), csaw:::.weightedFDR(pvals, rep(1, length(pvals)))))
pvals <- runif(1000)
stopifnot(all.equal(p.adjust(pvals, method="BH"), csaw:::.weightedFDR(pvals, rep(1, length(pvals)))))
pvals <- rep(1, 1000)
stopifnot(all.equal(p.adjust(pvals, method="BH"), csaw:::.weightedFDR(pvals, rep(1, length(pvals)))))
csaw:::.weightedFDR(numeric(0), numeric(0))

###################################################
# Reference results for the consolidation function; also tests clusterWindows and controlClusterFDR.
# There's not much point doing exact checks, because we'd just be re-implementing most of it.
# Note that there's no guaratee that the results and the reported FDR are the same, as this depends
# on the random p-values, so we don't check for equality.

checkResults <- function(data.list, result.list, pval.col=NULL, ..., true.pos) {
    out <- consolidateClusters(data.list, result.list, pval.col=pval.col, ...)

    # Checking that the clustering is fine.
    all.ids <- unlist(out$id)
    ref <- splitAsList(do.call(c, data.list), all.ids)
    stopifnot(all(unlist(range(ref))==out$region))

    # Checking that the right windows were chosen.
    if (is.null(pval.col)) { pval.col <- "PValue" }
    all.ps <- unlist(sapply(result.list, FUN=function(x) { x[,pval.col] }))
    was.sig <- !is.na(all.ids)
    if (any(was.sig) && any(!was.sig)) { stopifnot(max(all.ps[was.sig]) < min(all.ps[!was.sig])) }

    # Reporting the observed and estimated FDRs.
    np <- out$region[!overlapsAny(out$region, true.pos),]
    return(data.frame(Observed=length(np)/length(out$region), Estimated=out$FDR))
}

set.seed(100)
windows <- GRanges("chrA", IRanges(1:1000, 1:1000))
test.p <- runif(1000)
test.p[rep(1:2, 100) + rep(0:99, each=2) * 10] <- 0 

true.pos <- windows[test.p==0]
checkResults(list(windows), list(data.frame(PValue=test.p)), tol=0, target=0.05, true.pos=true.pos)
checkResults(list(windows), list(data.frame(PValue=test.p)), tol=10, target=0.05, true.pos=true.pos)

checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p), data.frame(PValue=test.p[1:10])), tol=0, target=0.05, true.pos=true.pos) # Multiple entries
checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p), data.frame(PValue=test.p[1:10])), equiweight=FALSE, tol=0, target=0.05, true.pos=true.pos)

# Smaller number of regions
set.seed(50)
test.p <- runif(1000)
test.p[rep(1:2, 50) + rep(0:49, each=2) * 10] <- 0  

true.pos <- windows[test.p==0]
checkResults(list(windows), list(data.frame(PValue=test.p)), tol=0, target=0.05, true.pos=true.pos)
checkResults(list(windows), list(data.frame(PValue=test.p)), tol=5, target=0.05, true.pos=true.pos)
checkResults(list(windows), list(data.frame(PValue=test.p)), tol=5, target=0.1, true.pos=true.pos)
checkResults(list(windows), list(data.frame(whee=test.p)), tol=2, pval.col="whee", target=0.05, true.pos=true.pos)

signs <- ifelse(rbinom(1000, 1, 0.5)!=0L, 1, -1)
checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p, logFC=signs), data.frame(PValue=test.p[1:10], logFC=signs[1:10])), 
             tol=0, target=0.05, true.pos=true.pos)
checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p, logFC=signs), data.frame(PValue=test.p[1:10], logFC=signs[1:10])), 
             tol=0, fc.col="logFC", target=0.05, true.pos=true.pos)

checkResults(list(windows), list(data.frame(PValue=test.p)), tol=0, grid.param=list(scale=5, iter=10), target=0.05, true.pos=true.pos) # Fiddling with grid search parameters.
checkResults(list(windows), list(data.frame(PValue=test.p)), tol=0, grid.param=list(len=11, it=10), target=0.05, true.pos=true.pos)

###################################################
