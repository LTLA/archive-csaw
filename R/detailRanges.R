detailRanges <- function(incoming, txdb, orgdb, dist=5000, promoter=c(3000, 0), max.intron=1e6)
# This gives three character vectors for each 'incoming'. The first specifies
# which features are wholly or partially overlapped by the current range.
# The second specifies what features are within 'tol' of the 5' end of 
# the incoming region, and the strand of that feature. The last specifies
# what features are within 'tol' of the 3' end of the incoming region.
# Promoters and exons are included. If there are multiple overlaps, the 
# first and last exon matching is shown.
#
# written by Aaron Lun
# 23 November, 2013
{
	require(GenomicFeatures)
	left.flank <- suppressWarnings(trim(flank(incoming, dist)))
	right.flank <- suppressWarnings(trim(flank(incoming, dist, start=FALSE)))

	# Obtain exons, and cleaning out the annotation.
	curex <- exonsBy(txdb, by="gene")
	curex <- GenomicRanges::unlist(curex)
	gene.id <- names(curex)
	gene.str <- as.logical(strand(curex)=="+")
	curex$exon_id <- NULL
	curex$exon_name <- NULL

	# Getting name annotation.
	anno <- select(orgdb, keys=gene.id, columns=c("SYMBOL"), keytype="ENTREZID")
	gene.name <- ifelse(is.na(anno$SYMBOL), paste0("ID:", anno$ENTREZID), anno$SYMBOL)
	
	# Splitting IDs, to avoid problems when genes are assigned to multiple locations.
	# The start uses '>' as we should be on the same gene by that stage; and exonBy
	# should give sorted locations w.r.t. start, if everything else is the same.
	gene.id <- as.integer(gene.id)
	is.diff <- c(TRUE, diff(gene.id)!=0L | diff(as.integer(seqnames(curex)))!=0L
			| diff(gene.str)!=0L | diff(start(curex)) > max.intron)
	gene.id <- cumsum(is.diff)
	ngenes <- sum(is.diff)

	# Assembling exon counts. Note that everything is sorted within each gene 
	# by exon start. We resort by exon end for the reverse strand to ensure
	# that '1' is the first exon in pathological cases. 
	output <- .Call("R_collate_exon_data", gene.id, gene.str, start(curex), end(curex), PACKAGE="csaw")
	if (is.character(output)) { stop(output) }
	ex.num <- output[[1]]

	##############################
	# Adding the gene bodies, using the extremes for each gene (returned as above).
	gb.collected <- output[[2]]
	gb.ref <- gb.collected[[1]]
	gb.ranges <- GRanges(seqnames(curex)[gb.ref], IRanges(gb.collected[[2]], gb.collected[[3]]),
		strand=!gene.str[gb.ref], seqinfo=seqinfo(curex))

	# Adding the promoter regions, based on the start/end of the gene bodies.
	if (length(promoter)!=2L) {  stop("need upstream/downstream specification in promoter argument") }
	prom.ranges <- suppressWarnings(promoters(gb.ranges, upstream=promoter[1], downstream=promoter[2]))

	# Expanding everything to account for the gene body and promoters.
	curex <- c(curex, prom.ranges, gb.ranges)
	gene.name <- c(gene.name, gene.name[gb.ref], gene.name[gb.ref])
	gene.id <- c(gene.id, gene.id[gb.ref], gene.id[gb.ref])
	ex.num <- c(ex.num, integer(length(gb.ref)), rep(-1L, length(gb.ref)))
	gene.str <- c(gene.str, gene.str[gb.ref], gene.str[gb.ref])
	
	###############################
	# Computing overlaps to everyone and his dog. Note that we don't do
	# intronic or promoter overlaps to the flanks; we don't care that
	# we're so-and-so base pairs away from the end of the promoter.
	full.lap <- findOverlaps(incoming, curex)
	flank.only <- ex.num > 0L
	left.lap <- findOverlaps(left.flank, curex[flank.only])
	left.dist <- start(incoming)[queryHits(left.lap)] - end(curex[flank.only])[subjectHits(left.lap)]
	right.lap <- findOverlaps(right.flank, curex[flank.only])
	right.dist <- start(curex[flank.only])[subjectHits(right.lap)] - end(incoming)[queryHits(right.lap)]  
	
	# Collating the left-overs.
	all.strs <- .Call("R_annotate", length(incoming), queryHits(full.lap)-1L, subjectHits(full.lap)-1L,
			queryHits(left.lap)-1L, which(flank.only)[subjectHits(left.lap)]-1L, left.dist,
			queryHits(right.lap)-1L, which(flank.only)[subjectHits(right.lap)]-1L, right.dist,
			gene.name, ex.num, gene.id, gene.str, PACKAGE="csaw")
	if (is.character(all.strs)) { stop(all.strs) }
	names(all.strs) <- c("overlap", "left", "right")
	return(all.strs)
}

