detailRanges <- function(incoming, txdb, orgdb, dist=5000, promoter=c(3000, 1000), 
    max.intron=1e6, key.field="ENTREZID", name.field="SYMBOL", ignore.strand=TRUE)
# This gives three character vectors for each 'incoming'. The first specifies
# which features are wholly or partially overlapped by the current range.
# The second specifies what features are within 'tol' of the 5' end of 
# the incoming region, and the strand of that feature. The last specifies
# what features are within 'tol' of the 3' end of the incoming region.
# Promoters and exons are included. If there are multiple overlaps, the 
# first and last exon matching is shown.
#
# written by Aaron Lun
# created 23 November 2013
# last modified 22 July 2015
{
	# Obtain exons, and cleaning out the annotation.
	curex <- exonsBy(txdb, by="gene")
	curex <- GenomicRanges::unlist(curex)
	gene.id <- names(curex)
	gene.str <- as.logical(strand(curex)=="+")
	curex$exon_id <- NULL
	curex$exon_name <- NULL

	# Getting name annotation.
	summarized <- rle(gene.id)
	anno <- suppressMessages(select(orgdb, keys=summarized$values, columns=name.field, keytype=key.field))
	re.summarized <- rle(as.character(anno[[key.field]]))
	if (any(re.summarized$length!=1L)) { 
		warning("many-to-one relationship between key and name fields, using first value only") 
		first.of.each <- c(1L, cumsum(re.summarized$length)[-length(re.summarized$length)] + 1L)
		anno <- anno[first.of.each,]
	}

	all.names <- list()	
	do.check <- !key.field %in% name.field # Redundant, if it's already being reported.
	for (x in seq_along(name.field)) {
		cur.name <- inverse.rle(list(values=anno[[name.field[x]]], lengths=summarized$length))
		if (!is.character(cur.name)) { cur.name <- as.character(cur.name) } 
		if (do.check) { 
			failed <- is.na(cur.name)
			if (any(failed)) { cur.name[failed] <- paste0("<", gene.id[failed], ">") }
		}
		all.names[[x]] <- cur.name
	}
	gene.name <- do.call(paste, c(all.names, sep=";"))
	
	# Splitting IDs, to avoid problems when genes are assigned to multiple locations.
	# exonBy should give sorted locations by chr/strand/start/end, so gap between 
	# start of each exon and the end of the previous one should give the intron length (ignore nesting).
	nexons <- length(gene.id)
	is.diff <- c(TRUE, gene.id[-1]!=gene.id[-nexons] | diff(as.integer(seqnames(curex)))!=0L
		| diff(gene.str)!=0L | start(curex)[-1] - end(curex)[-nexons] - 1L > max.intron)
	gene.id <- cumsum(is.diff)
	ngenes <- sum(is.diff)

	# Assembling exon counts. Note that everything is sorted within each gene 
	# by exon start. We resort by exon end for the reverse strand to ensure
	# that '1' is the first exon in pathological cases. 
	output <- .Call(cxx_collate_exon_data, gene.id, gene.str, start(curex), end(curex))
	if (is.character(output)) { stop(output) }
	ex.num <- output[[1]]

	##############################
	# Adding the gene bodies, using the extremes for each gene (returned as above).
	gb.collected <- output[[2]]
	gb.ref <- gb.collected[[1]]
	gb.ranges <- GRanges(seqnames(curex)[gb.ref], IRanges(gb.collected[[2]], gb.collected[[3]]),
		strand=!gene.str[gb.ref], seqinfo=seqinfo(curex))
	names(gb.ranges) <- names(curex)[gb.ref]

	# Adding the promoter regions, based on the start/end of the gene bodies.
	if (length(promoter)!=2L) {  stop("need upstream/downstream specification in promoter argument") }
	prom.ranges <- suppressWarnings(trim(promoters(gb.ranges, upstream=promoter[1], downstream=promoter[2])))
	names(prom.ranges) <- names(gb.ranges)

	# Expanding everything to account for the gene body and promoters.
	curex <- c(curex, prom.ranges, gb.ranges)
	gene.name <- c(gene.name, gene.name[gb.ref], gene.name[gb.ref])
	gene.id <- c(gene.id, gene.id[gb.ref], gene.id[gb.ref])
	ex.num <- c(ex.num, integer(length(gb.ref)), rep(-1L, length(gb.ref)))
	gene.str <- c(gene.str, gene.str[gb.ref], gene.str[gb.ref])
	
	# Returning the useful stuff, if no overlaps are requested.
	if (missing(incoming)) { 
		curex$symbol <- gene.name
		curex$exon <- ex.num
		curex$internal <- gene.id
		return(curex)
	}

	###############################
	# Computing overlaps to everyone and his dog. Note that we don't do
	# intronic or promoter overlaps to the flanks; we don't care that
	# we're so-and-so base pairs away from the end of the promoter. We also
	# rule out negative or zero distances i.e. those that would overlap the
	# region itself (zero distance means 1-based end and start are equal).
	# We may or may not consider strandedness (flanking sets ignore.true
	# so that we always get left/right flanks).

	full.lap <- findOverlaps(incoming, curex, ignore.strand=ignore.strand)
	flank.only <- ex.num > 0L
	to.flank <- curex[flank.only]

	left.flank <- suppressWarnings(trim(flank(incoming, dist, ignore.strand=TRUE)))
	left.lap <- findOverlaps(left.flank, to.flank, ignore.strand=ignore.strand)
	left.dist <- start(incoming)[queryHits(left.lap)] - end(to.flank)[subjectHits(left.lap)]
	left.nolap <- left.dist > 0L
	left.lap <- left.lap[left.nolap,]
	left.dist <- left.dist[left.nolap]

	right.flank <- suppressWarnings(trim(flank(incoming, dist, start=FALSE, ignore.strand=TRUE)))
	right.lap <- findOverlaps(right.flank, to.flank, ignore.strand=ignore.strand)
	right.dist <- start(to.flank)[subjectHits(right.lap)] - end(incoming)[queryHits(right.lap)]  
	right.nolap <- right.dist > 0L
	right.lap <- right.lap[right.nolap,]
	right.dist <- right.dist[right.nolap]
	
	# Collating the left-overs.
	all.strs <- .Call(cxx_annotate_overlaps, length(incoming), queryHits(full.lap)-1L, subjectHits(full.lap)-1L,
			queryHits(left.lap)-1L, which(flank.only)[subjectHits(left.lap)]-1L, left.dist,
			queryHits(right.lap)-1L, which(flank.only)[subjectHits(right.lap)]-1L, right.dist,
			gene.name, ex.num, gene.id, gene.str)
	if (is.character(all.strs)) { stop(all.strs) }
	names(all.strs) <- c("overlap", "left", "right")
	return(all.strs)
}

