# This tests that the empirical FDR performs its calculations correctly.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))

comp <- function(total.n, n.false, n.true, weights=NULL) {
    n.clusters <- n.false + n.true
    merged.ids <- sample(n.clusters, total.n, replace=TRUE)
    tab <- data.frame(logFC=abs(rnorm(total.n)), logCPM=runif(total.n, -2, 1),
                      PValue=rbeta(total.n, 1, 100))
    
    # Adding true nulls.
    is.false <- merged.ids <= n.false
    tab$PValue[is.false] <- runif(sum(is.false))
    tab$logFC[is.false] <- rnorm(sum(is.false))

    out <- empiricalFDR(merged.ids, tab, weight=weights)
    stopifnot(identical(rownames(out), as.character(sort(unique(merged.ids)))))

    # Checking calculations:    
    new.p <- tab$PValue/2
    new.p[tab$logFC < 0] <- 1 - new.p[tab$logFC < 0]
    tab2 <- tab
    tab2$PValue <- new.p
    ref <- combineTests(merged.ids, tab2, weight=weights)
    stopifnot(all(abs(out$PValue - ref$PValue) <= 1e-6))

    alt <- empiricalFDR(merged.ids, tab, weight=weights, neg.down=FALSE)
    tab2 <- tab
    tab2$PValue <- 1-new.p
    ref <- combineTests(merged.ids, tab2, weight=weights)
    stopifnot(all(abs(alt$PValue - ref$PValue) <= 1e-6))

    emp.fdr <- findInterval(out$PValue, sort(alt$PValue))/rank(out$PValue, ties.method="max")
    emp.fdr <- pmin(emp.fdr, 1)
    o <- order(out$PValue, decreasing=TRUE)
    emp.fdr[o] <- cummin(emp.fdr[o])
    stopifnot(all(abs(emp.fdr - out$FDR) <= 1e-6))
    
    # Calculating the actual FDR.
    is.sig <- out$FDR <= 0.05
    is.true <- as.integer(rownames(out)) > n.false
    actual.fdr <- sum(is.sig & !is.true)/sum(is.sig)
    return(actual.fdr)
}

set.seed(1000)

# All of these values are pretty close to 5%, though it does get a bit unstable with lower numbers of tests.
comp(10000, 5000, 1000)
comp(10000, 5000, 1000)
comp(10000, 5000, 1000)
comp(10000, 5000, 1000)

comp(5000, 500, 100)
comp(5000, 500, 100)
comp(5000, 500, 100)
comp(5000, 500, 100)

comp(1000, 500, 100)
comp(1000, 500, 100)
comp(1000, 500, 100)
comp(1000, 500, 100)

# Repeating with weights.
comp(10000, 5000, 1000, weights=runif(10000))
comp(5000, 500, 100, weights=runif(5000))
comp(1000, 500, 100, weights=runif(1000))

###################################################################################################
# Checking for sane behaviour when empty.

empiricalFDR(integer(0), data.frame(PValue=numeric(0), logCPM=numeric(0), logFC=numeric(0)), weight=numeric(0))

###################################################################################################
# End.


