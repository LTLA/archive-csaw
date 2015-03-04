# This tests the behaviour of the detailRanges function, by comparing the output of the
# strings to the output of the function proper.

suppressPackageStartupMessages(require(csaw))
suppressPackageStartupMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene))
suppressPackageStartupMessages(require(org.Mm.eg.db))
source("simsam.R")
require(parallel)

########################################################################################
# Checking the sensibility of the exon numbering, the promoters and gene bodies, in each case.

checkranges <- function(ref, up, down) {
	by.gene <- split(ref, ref$internal)
	mclapply(names(by.gene), FUN=function(x) { 
		entry <- by.gene[[x]]
		all.exons <- entry[entry$exon >= 1L]
		true.gb <- range(all.exons)
		cur.gb <- entry[entry$exon==-1L]
		names(cur.gb) <- elementMetadata(cur.gb) <- NULL
		if (!identical(cur.gb, true.gb)) { 
			print(all.exons)
			print(cur.gb)
			print(true.gb)
			stop("gene body identification is wrong") }

		promoter <- entry[entry$exon==0L]
		first.exon <- entry[entry$exon==1L]
		if (!identical(seqnames(promoter), seqnames(first.exon)) || !identical(strand(promoter), strand(first.exon))) {
			stop("promoter characteristics are not consistent with exons") }
		chrlen <- seqlengths(ref)[[runValue(seqnames(promoter))]]

		if (as.character(runValue(strand(entry)))=="+") { 
			if (is.unsorted(start(all.exons)[order(all.exons$exon)])) { 
				print(all.exons)
				stop("exon numbering is wrong for forward-strand gene") 
			}
			if (max(1L, start(first.exon) - up) != start(promoter) || min(chrlen, start(first.exon) + down - 1L) != end(promoter)) {
				print(promoter)
				print(first.exon)
				stop("promoter generation is wrong for forward-strand gene") }
		} else {
			if (is.unsorted(end(all.exons)[order(all.exons$exon, decreasing=TRUE)])) { 
				print(all.exons)
				stop("exon numbering is wrong for reverse strand") 
			}
			if (max(1L, end(first.exon) - down + 1L) != start(promoter) || min(chrlen, end(first.exon) + up) != end(promoter)) {
				print(promoter)
				print(first.exon)
				stop("promoter generation is wrong for reverse-strand gene") }
		}
		return(0)
	}, mc.cores=8)
	return(ref)
}

up <- 3000
down <- 1000
ref <- detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, promoter=c(up, down))
checkranges(ref, up, down)

########################################################################################
### Making a comparator function, to check proper string construction.

comp <- function(incoming, ref, dist=5000) {
	olap <- findOverlaps(incoming, ref)	
	anno <- detailRanges(incoming, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, dist=dist)	

	.getmode <- function(collected.modes) {
		collected.modes <- sort(collected.modes)
		if (collected.modes[1]==-1L) {
			if (length(collected.modes)==1L) { curmode <- "I" } 
			collected.modes <- collected.modes[-1]
		}
		if (length(collected.modes)) {
			consec.start <- which(c(TRUE, diff(collected.modes)!=1L))
			consec.end <- c(consec.start[-1] - 1L, length(collected.modes))
			all.strings <- ifelse(consec.start==consec.end, collected.modes[consec.start], paste0(collected.modes[consec.start], "-", collected.modes[consec.end]	))
			curmode <- paste(all.strings, collapse=",")
		}
		return(curmode)
	}

	# Checking overlaps.
	relevants <- split(subjectHits(olap), queryHits(olap))
	test.anno <- mclapply(names(relevants), FUN=function(it) { 
		actual.index <- as.integer(it)
		by.gene <- split(relevants[[it]], ref$internal[relevants[[it]]])
		cur.anno <- anno$overlap[actual.index]

		# Assembling the string based on what's going on.
		for (g in names(by.gene)) { 
			collected.modes <- ref$exon[by.gene[[g]]]
			curmode <- .getmode(collected.modes)
			fstring <- paste(ref$symbol[by.gene[[g]][1]], curmode, strand(ref[by.gene[[g]][1]]), sep="|") 
			if (!grepl(fstring, cur.anno, fixed=TRUE)) {
				print(fstring)
				print(cur.anno)
				stop("could not find overlap")
			} 
			cur.anno <- sub(fstring, "", cur.anno, fixed=TRUE)
		}
		return(cur.anno)
	}, mc.cores=8)
	stopifnot(all(nchar(gsub(",", "", unlist(test.anno)))==0L))

	# Checking left and right overlaps.
	for (mode in 1:2) { 
		if (mode==1L) { 
			test.anno <- anno$left
			olap <- findOverlaps(GRanges(seqnames(incoming), IRanges(start(incoming)-dist, start(incoming)-1L)), ref)	
			relevant.x <- split(subjectHits(olap), queryHits(olap)) 
		} else {
			test.anno <- anno$right
			olap <- findOverlaps(GRanges(seqnames(incoming), IRanges(end(incoming)+1L, end(incoming)+dist)), ref)	
			relevant.x <- split(subjectHits(olap), queryHits(olap)) 
		}

		new.anno <- mclapply(names(relevant.x), FUN=function(it) { 
			actual.index <- as.integer(it)
			relevant.x[[it]] <- setdiff(relevant.x[[it]], relevants[[it]]) # Getting rid of the centers.
			by.gene <- split(relevant.x[[it]], ref$internal[relevant.x[[it]]])
			cur.anno <- test.anno[actual.index]

			# Assembling the string based on what's going on.
			for (g in names(by.gene)) { 
				chosen <- ref[by.gene[[g]]]
				collected.modes <- chosen$exon
				keep <- collected.modes > 0L
				collected.modes <- collected.modes[keep]
				chosen <- chosen[keep]
				if (!any(keep)) { next }

				curmode <- .getmode(collected.modes)
				fstring <- paste(chosen$symbol[1], curmode, strand(chosen[1]), sep="|")
				all.dist <- pmax(start(incoming[actual.index]), start(chosen)) - pmin(end(incoming[actual.index]), end(chosen))
				if (any(all.dist <= 0L || all.dist > dist)) { stop("distances out of range") }
				fstring <- paste0(fstring, "[", min(all.dist), "]")
					
				if (!grepl(fstring, cur.anno, fixed=TRUE)) {
					print(fstring)
					print(cur.anno)
					stop("could not find flank")
				} 
				cur.anno <- sub(fstring, "", cur.anno, fixed=TRUE)
			}
			return(cur.anno)
		}, mc.cores=8)
		stopifnot(all(nchar(gsub(",", "", unlist(new.anno)))==0L))
	}
	
	return(data.frame(O=sum(nchar(anno$overlap)!=0L), L=sum(nchar(anno$left)!=0L), R=sum(nchar(anno$right)!=0L)))
}

set.seed(1847382)
chromos <- seqlengths(ref)
chromos <- chromos[chromos > 1e7]

all.win <- generateWindows(chromos*1.8, 1e3, 1000)
comp(all.win, ref)

all.win <- generateWindows(chromos*1.8, 1e3, 2000)
comp(all.win, ref, dist=2000)

all.win <- generateWindows(chromos*1.8, 1e3, 2000)
comp(all.win, ref, dist=10000)

########################################################################################
#### Checking key, name field options.

output <- detailRanges(all.win, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
		    orgdb=org.Mm.eg.db, name.field=c("ENTREZID"))
head(output$overlap, 30)

output <- detailRanges(all.win, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
	    orgdb=org.Mm.eg.db, name.field=c("SYMBOL", "ENTREZID"))
head(output$overlap, 30)

suppressPackageStartupMessages(require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene))
suppressPackageStartupMessages(require(org.Sc.sgd.db))

allr <-detailRanges(txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, orgdb=org.Sc.sgd.db, key.field='ORF', name.field='GENENAME')
allr

########################################################################################
########################################################################################
########################################################################################
# End.
