.normalizeChIP <- function(counts, lib.sizes=colSums(counts), refColumn=NULL, prior.count=0.1, logratioTrim = 0.3, sumTrim = 0.05)
# This performs TMM normalization with some modifications to deal with low 
# counts and other problems specific to ChIP-seq data. It also has a slight
# difference in implementation to improve behaviour around ties.
#
# written by Aaron Lun
# 7 October 2013
{ 
   	# Adding a prior count to bring (0, n) pairs into the trimming framework.
	if (prior.count <= 0) { stop("prior count should be positive") }
	expr<-cpm(counts, lib.size=lib.sizes, log=TRUE, prior.count=prior.count)

	# Choosing the reference based on the most "middle" column (i.e. easiest to compare to all others).
	if (is.null(refColumn)) { 
		f<-apply(expr, 2, FUN=median)
		refColumn<-which.min(abs(f-mean(f)))
	}
	if (sumTrim < 0 || sumTrim > 0.5 || logratioTrim < 0 || logratioTrim > 0.5) { 
		stop("trimming proportions should lie between 0 and 0.5");
	}

	# Iterating across and using the trimmed mean of M-values.
	nlibs<-ncol(counts)
	normfacs<-numeric(nlibs)
	for (lib in 1:nlibs) {
		if (lib==refColumn) { next }
	
		# Discarding (0, 0) pairs with no information for a pairwise comparison. 
		keep<-counts[,lib]!=0L & counts[,refColumn]!=0L
		o<-expr[keep,lib]
		r<-expr[keep,refColumn]
		n<-sum(keep)
		
		# Trimming off extreme abundances with 'order' instead of 'rank', avoid issues with ties.
		lo <- floor(n * sumTrim) + 1L
		hi <- n + 1L - lo
		a<-r+o
		akeep<-logical(n)
		akeep[order(a)[lo:hi]]<-TRUE
		
		# Trimming of extreme M-values.
		lo <- floor(n * logratioTrim) + 1L
		hi <- n + 1L - lo
		m<-o-r
		mkeep<-logical(n)
		mkeep[order(m)[lo:hi]]<-TRUE

		# Computing the trimmed mean of M-values. Precision weighting is not performed to avoid 
		# upweighting high abundance bins which are more likely to be (differentially) enriched.
		if (!any(mkeep & akeep)) { stop("trimming proportions are too aggressive") }
		normfacs[lib]<-mean(m[mkeep & akeep])	
	}

	# Centering the normalization factors around a value of 1.
	normfacs<-normfacs-mean(normfacs)
	return(2^normfacs)
}
 
normalizeChIP <- function(counts, lib.sizes=colSums(counts), weighted=FALSE, ...) 
# This is a wrapper to perform TMM normalization with non-standard library 
# sizes (e.g. due to filtering) and weighting turned off. Purely for 
# convenience, as I don't want to put in an extra three lines per call 
# (especially when there are multiple calls to test robustness).
#
# written by Aaron Lun
# 19 November, 2013
{
	y<-DGEList(counts, lib.size=lib.sizes)
	y<-calcNormFactors(y, doWeighting=weighted, ...)
	y$samples$norm.factors
}
