\name{annotateRange}
\alias{annotateRange}

\title{Add annotation to ranges}
\description{Add exon-based annotation to specified regions}

\usage{
annotateRange(incoming, txdb, orgdb, dist=5000, promoter=3000, max.intron=1e6)
}

\arguments{
\item{incoming}{a \code{GRanges} object containing the ranges to be annotated}
\item{txdb}{a TranscriptDb object for the genome of interest}
\item{orgdb}{a Genome wide annotation object for the genome of interest}
\item{dist}{an integer scalar specifying the flanking distance to annotate}
\item{promoter}{an integer scalar specifying the size of a promoter region}
\item{max.intron}{an integer scalar indicating the maximum distance between exons}
}

\value{
A list with \code{overlap}, \code{left} and \code{right} elements, each a character
vector of length equal to the number of ranges in \code{incoming}. Each non-empty
string records the gene symbol, the overlapped exons and the strand. For \code{left}
and \code{right}, the gap between the range and the annotated feature is also included.
}

\details{
For annotated features overlapping a region, the character string in the
\code{overlap} output vector will be \code{GENE|EXONS|STRAND}. The \code{GENE}
is the gene symbol, with an alternative as \code{ID:XXX} for the Entrez ID if
the symbol is unavailable. The \code{EXONS} indicate the exon or range of exons
that are overlapped. The \code{STRAND} is, obviously, the strand on which the
gene is coded. For annotated regions flanking the region within \code{dist},
the character string in the \code{left} or \code{right} output vectors will
have an additional \code{DIST} value. This represents the gap between the edge
of the region and the closest exon for that gene.

Exons are numbered in order of increasing start or end position for genes
on the forward or reverse strands, respectively. Promoters are defined as
the region of length \code{promoter} upstream of the gene TSS, itself 
defined as the start of the first exon (for genes on the forward strand)
or the end of the last exon (otherwise). All promoters are marked as
exon 0 for simplicity. Exon ranges in \code{EXON} are reported from
the lowest to highest numbering exons. If the region overlaps an intron, 
it is labelled with \code{I} in \code{EXON}. No overlaps to introns are
reported if there is an overlap to an exon.

Note that promoter and intronic annotations are only reported for the
\code{overlap} vector. This reduces redundancy in the output. For example, it
makes little sense to report that the region is both flanking and overlapping
an intron. Similarly, the value of \code{DIST} is also more relevant when it is
reported to the nearest exon rather than to an intron (in which case, the
distance would be zero if the intron overlaps the region). In cases where
the distance is reported to the first exon, it can be used to refine the
choice of \code{promoter}.

The \code{max.intron} value is necessary to deal with genes that have ambiguous
locations on the genome. If a gene has exons on different chromosomes, its
location is uncertain and the gene is partitioned into two sets of exons for
separate processing. However, this is less obvious when the ambiguous locations
belong to the same chromosome. The \code{max.intron} value protects against
excessively large genes that may occur from considering those locations as 
a single transcriptional unit.
}

\author{Aaron Lun}

\examples{ 
require(org.Mm.eg.db)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)

current <- readRDS(system.file("exdata", "exrange.rds", package="csaw"))
output <- annotateRange(current, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db)
head(output$overlap)
head(output$right)
head(output$left)
}
