# This tests the behaviour of the detailRanges function, by comparing the output of the
# strings to the output of the function proper.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))
suppressWarnings(suppressPackageStartupMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene)))
suppressWarnings(suppressPackageStartupMessages(require(org.Mm.eg.db)))
source("simsam.R")

########################################################################################
# Checking the sensibility of the exon numbering, the promoters and gene bodies, in each case.

checknames <- function(ref) {
	gene.ids <- names(ref)
	gene.sym <- ref$symbol
	nuniq1 <- sapply(split(gene.ids,gene.sym), FUN=function(x) { length(unique(x)) })
	nuniq2 <- sapply(split(gene.sym,gene.ids), FUN=function(x) { length(unique(x)) })
	if (any(nuniq1!=1L) || any(nuniq2!=1L)) { stop("name to symbol assignment is not unique") }
	return(invisible(NULL))
}

checkranges <- function(ref, up, down) {
	exonic <- ref[ref$exon>=1L]
	gb <- unlist(range(split(exonic, exonic$internal)))
	test.gb <- ref[ref$exon==-1L]
	names(test.gb) <- test.gb$internal
	elementMetadata(test.gb) <- NULL
	if (!identical(gb, test.gb)) { stop("differences in gene body identification") }

	promoters <- ref[ref$exon==0L]
 	first.exon <- ref[ref$exon==1L]
	is.forward <- as.logical(strand(first.exon)=="+")
	new.prom <- first.exon
	suppressWarnings(start(new.prom) <- ifelse(is.forward, start(first.exon) - up, end(first.exon) - down + 1L))
	suppressWarnings(end(new.prom) <- ifelse(is.forward, start(first.exon) + down - 1L, end(first.exon) + up))
	new.prom <- trim(new.prom)
	new.prom$exon <- 0L
	if (!identical(new.prom, promoters)) { stop("differences in promoter identification") }

	forward.exons <- exonic[strand(exonic)=="+"]
	o <- order(forward.exons$internal, forward.exons$exon)
	n <- length(o)
	forward.exons <- forward.exons[o]
	out.of.order <- c(FALSE, forward.exons$internal[-1]==forward.exons$internal[-n] & 
		start(forward.exons)[-n] > start(forward.exons)[-1])
	if (any(out.of.order)) { stop("exon ranking for forward-strand genes is incorrect") }

	reverse.exons <- exonic[strand(exonic)=="-"]
	o <- order(reverse.exons$internal, reverse.exons$exon)
	n <- length(o)
	reverse.exons <- reverse.exons[o]
	out.of.order <- c(FALSE, reverse.exons$internal[-1]==reverse.exons$internal[-n] & 
		end(reverse.exons)[-n] < end(reverse.exons)[-1])
	if (any(out.of.order)) { stop("exon ranking for reverse-strand genes is incorrect") }

	checknames(ref)	
	return(promoters)
}

up <- 3000
down <- 1000
ref <- detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, promoter=c(up, down))
checkranges(ref, up, down)

up <- 2000
down <- 500
ref <- detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, promoter=c(up, down))
checkranges(ref, up, down)

up <- 5000
down <- 0
ref <- detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, promoter=c(up, down))
checkranges(ref, up, down)

########################################################################################
### Setting up a comparator function, in C++.

library(inline)

# Assumes sorting by 'queries' and then 'subjects', zero-indexing for both.
# This is basically a rehash of the C++ code in the csaw package; running
# it in R is unacceptably slow as we need loops within loops (no easy vectorization
# to construct each string). My hope is that two separate implementations can't be 
# wrong, or won't stuff up at the same time.

compiled <- cxxfunction(sig=c(id="integer", name="character", exon="integer", strand="character", 
		nquery="integer", queries="integer", subjects="integer", distances="integer"), body='
const int* gids=INTEGER(id);
const int* exonnum=INTEGER(exon);
const int* qid=INTEGER(queries);
const int nqid=LENGTH(queries);
const int* sid=INTEGER(subjects);
const int* dists=INTEGER(distances);
const bool dodist=LENGTH(distances)>0;

SEXP output=PROTECT(allocVector(STRSXP, INTEGER(nquery)[0]));
for (int i=0; i<nqid; ++i) { SET_STRING_ELT(output, i, mkChar("")); }

int counter=0, lastpos=0;
while (counter < nqid) { 
	lastpos=counter;
	std::map<int, std::deque<std::pair<int, int> > > collected;
    do { 
		collected[gids[sid[counter]]].push_back(std::make_pair(exonnum[sid[counter]], counter));
		++counter; 
	} while (counter < nqid && qid[counter-1]==qid[counter]);
	
	// Manufacturing a string.
	std::stringstream out;
	for (std::map<int, std::deque<std::pair<int, int> > >::iterator itc=collected.begin(); itc!=collected.end(); ++itc) {
		int holding=-1;
		std::deque<std::pair<int, int> >& current=itc->second;
		std::sort(current.begin(), current.end());
		out << CHAR(STRING_ELT(name, sid[current.front().second])) << "|";

		// Deciding whether to add an intron, or skip to the next non-intronic element.
		if (current.front().first==-1) {
			if (current.size()==1) {  out << "I"; } 
			else {
				current.pop_front(); 
				out << current.front().first;
			} 
		} else { out << current.front().first; }
			
		// Constructing the full string.
		for (size_t c=1; c<current.size(); ++c) { 
			if (current[c-1].first+1==current[c].first) { holding=current[c].first; } 
			else {
				if (holding >= 0) {
					out << "-" << holding;
					holding=-1;
				}
				out << "," << current[c].first;
			}
		}

		if (holding >= 0) { out << "-" << holding; }
		out << "|" << CHAR(STRING_ELT(strand, sid[current.front().second]));
			
		// Taking the distance to the first element (assuming introns and promoters have been tossed out).
		if (dodist) { 
			std::deque<int> temp(current.size());
			for (size_t c=0; c<current.size(); ++c) { temp[c]=dists[current[c].second]; }
			out << "[" << *std::min_element(temp.begin(), temp.end()) << "]"; 
		}
		out << ",";  
	}

	// Getting rid of the loose comma.
	std::string blah=out.str();
	blah.erase(blah.size()-1);
	SET_STRING_ELT(output, qid[counter-1], mkChar(blah.c_str()));
}

UNPROTECT(1);
return output;
', includes="#include <sstream>\n#include <deque>\n#include <map>\n#include <algorithm>")

########################################################################################
### Making a comparator function, to check proper string construction.

comp <- function(incoming, up, down, dist=5000, ignore.strand=TRUE) {
	ref <- detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, promoter=c(up, down))
	olap <- findOverlaps(incoming, ref)	
	anno <- detailRanges(incoming, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, dist=dist, 
		promoter=c(up, down), ignore.strand=ignore.strand)

	# Checking overlaps.
	olap <- findOverlaps(incoming, ref, ignore.strand=ignore.strand)	
	obs <- compiled(ref$internal, ref$symbol, ref$exon, as.character(strand(ref)),
			length(incoming), queryHits(olap)-1L, subjectHits(olap)-1L, integer(0))
	stopifnot(identical(anno$overlap, obs))

	# Checking flanks.
	for (mode in 1:2) { 
		if (mode==1L) { 
			test.anno <- anno$left
			olap <- findOverlaps(GRanges(seqnames(incoming), IRanges(start(incoming)-dist, start(incoming)-1L), 
				strand=strand(incoming)), ref, ignore.strand=ignore.strand)
			relative.dist <- start(incoming)[queryHits(olap)] - end(ref)[subjectHits(olap)]
		} else {
			test.anno <- anno$right
			olap <- findOverlaps(GRanges(seqnames(incoming), IRanges(end(incoming)+1L, end(incoming)+dist),
				strand=strand(incoming)), ref, ignore.strand=ignore.strand)
			relative.dist <- start(ref)[subjectHits(olap)] - end(incoming)[queryHits(olap)]
		}

		keep <- relative.dist > 0L &  ref$exon[subjectHits(olap)] > 0L
		olap <- olap[keep,]
		relative.dist <- relative.dist[keep]
		obs <- compiled(ref$internal, ref$symbol, ref$exon, as.character(strand(ref)),
				length(incoming), queryHits(olap)-1L, subjectHits(olap)-1L, relative.dist)
		stopifnot(identical(obs, test.anno))
	}
	return(c(O=sum(nchar(anno$overlap)!=0L), L=sum(nchar(anno$left)!=0L), R=sum(nchar(anno$right)!=0L)))
}

set.seed(1847382)
chromos <- seqlengths(ref)
chromos <- chromos[chromos > 1e7]

all.win <- generateWindows(chromos*1.8, 1e2, 1000)
comp(all.win, up=3000, down=1000)

all.win <- generateWindows(chromos*1.8, 1e2, 2000)
comp(all.win, up=3000, down=1000, dist=2000)

all.win <- generateWindows(chromos*1.8, 1e2, 2000)
comp(all.win, up=5000, down=500, dist=10000)

# Checking behaviour upon strandedness.
strand(all.win) <- sample(c("+", "-"), length(all.win), replace=TRUE) 
comp(all.win, up=5000, down=500, dist=10000, ignore.strand=FALSE)
comp(all.win, up=2000, down=1000, dist=5000, ignore.strand=FALSE)
strand(all.win) <- "*"

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

allr <- detailRanges(txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, orgdb=org.Sc.sgd.db, key.field='ORF', name.field='GENENAME')
checknames(allr)
allr

allr <- detailRanges(txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, orgdb=org.Sc.sgd.db, key.field='ORF', name.field=c('GENENAME', 'ORF'))
checknames(allr)
allr

########################################################################################
########################################################################################
########################################################################################
# End.
