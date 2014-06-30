\name{detailRanges}
\alias{detailRanges}

\title{Add annotation to ranges}
\description{Add detailed exon-based annotation to specified genomic regions.}

\usage{
detailRanges(incoming, txdb, orgdb, dist=5000, 
    promoter=c(3000, 1000), max.intron=1e6)
}

\arguments{
\item{incoming}{a \code{GRanges} object containing the ranges to be annotated}
\item{txdb}{a \code{TranscriptDb} object for the genome of interest}
\item{orgdb}{a genome wide annotation object for the genome of interest}
\item{dist}{an integer scalar specifying the flanking distance to annotate}
\item{promoter}{an integer vector of length 2 defining the promoter as some distance upstream and downstream from the TSS}
\item{max.intron}{an integer scalar indicating the maximum distance between exons}
}

\value{
If \code{incoming} is not provided, a \code{GRanges} object will be returned
containing ranges for the exons, promoters and gene bodies. Gene IDs, symbol
and exon numbers are also stored as metadata.

If \code{incoming} is a \code{GRanges} object, a list will be returned with
\code{overlap}, \code{left} and \code{right} elements. Each element is a
character vector of length equal to the number of ranges in \code{incoming}.
Each non-empty string records the gene symbol, the overlapped exons and the
strand. For \code{left} and \code{right}, the gap between the range and the
annotated feature is also included.
}

\details{
This function adds exon-based annotations to a given set of genomic regions, in
the form of compact character strings specifying the features overlapping and
flanking each region. The aim is to determine the genic context of empirically
identified regions. This allows some basic biological interpretation of
binding/marking in those regions. All neighbouring genes within a specified
range are reported, rather than just the closest gene to the region.

For annotated features overlapping a region, the character string in the
\code{overlap} output vector will be of the form \code{GENE|EXONS|STRAND}. The
\code{GENE} is the gene symbol, with an alternative as \code{ID:XXX} for the
Entrez ID if the symbol is unavailable. The \code{EXONS} indicate the exon or
range of exons that are overlapped. The \code{STRAND} is, obviously, the strand
on which the gene is coded. For annotated regions flanking the region within
\code{dist}, the character string in the \code{left} or \code{right} output
vectors will have an additional \code{DIST} value. This represents the gap
between the edge of the region and the closest exon for that gene.

Exons are numbered in order of increasing start or end position for genes
on the forward or reverse strands, respectively. Promoters are defined as
the region of length \code{promoter} upstream of the gene TSS, itself 
defined as the start of the first exon (for genes on the forward strand)
or the end of the last exon (otherwise). All promoters are marked as
exon 0 for simplicity. Exon ranges in \code{EXON} are reported from as a
comma-separated list where stretches of consecutive exons are summarized into a
range. If the region overlaps an intron, it is labelled with \code{I} in
\code{EXON}. No intronic overlaps are reported if there is an exonic overlap.

Note that promoter and intronic annotations are only reported for the
\code{overlap} vector to reduce redundancy in the output. For example, it makes
little sense to report that the region is both flanking and overlapping an
intron. Similarly, the value of \code{DIST} is more relevant when it is
reported to the nearest exon rather than to an intron (in which case, the
distance would be zero if the intron overlaps the region). In cases where the
distance is reported to the first exon, it can be used to refine the choice of
\code{promoter}.

The \code{max.intron} value is necessary to deal with genes that have ambiguous
locations on the genome. If a gene has exons on different chromosomes, its
location is uncertain and the gene is partitioned into two sets of exons for
separate processing. However, this is less obvious when the ambiguous locations
belong to the same chromosome. The \code{max.intron} value protects against
excessively large genes that may occur from considering those locations as 
a single transcriptional unit.

If \code{incoming} is missing, then the annotation will be provided directly to
the user in the form of a \code{GRanges} object. This may be more useful when
further work on the annotation is required. Exon numbers are provided in the 
metadata with promoters and gene bodies labelled as 0 and -1, respectively.
Overlaps to introns can be identified by finding those regions that overlap 
with gene bodies but not with any of the corresponding exons.
}

\author{Aaron Lun}

\examples{ 
require(org.Mm.eg.db)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)

current <- readRDS(system.file("exdata", "exrange.rds", package="csaw"))
output <- detailRanges(current, TxDb.Mmusculus.UCSC.mm10.knownGene, 
    org.Mm.eg.db)
head(output$overlap)
head(output$right)
head(output$left)

output <- detailRanges(current, TxDb.Mmusculus.UCSC.mm10.knownGene, 
    org.Mm.eg.db, promoter=c(2000, 1000))
head(output$overlap)
head(output$right)
head(output$left)

detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db)
}

\keyword{annotation}
