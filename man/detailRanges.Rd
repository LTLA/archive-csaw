\name{detailRanges}
\alias{detailRanges}

\title{Add annotation to ranges}
\description{Add detailed exon-based annotation to specified genomic regions.}

\usage{
detailRanges(incoming, txdb, orgdb, dist=5000, promoter=c(3000, 1000), 
    max.intron=1e6, key.field="ENTREZID", name.field="SYMBOL", 
    ignore.strand=TRUE)
}

\arguments{
\item{incoming}{a \code{GRanges} object containing the ranges to be annotated}
\item{txdb}{a \code{TranscriptDb} object for the genome of interest}
\item{orgdb}{a genome wide annotation object for the genome of interest}
\item{dist}{an integer scalar specifying the flanking distance to annotate}
\item{promoter}{an integer vector of length 2, where first and second values define the promoter as some distance upstream and downstream from the TSS, respectively}
\item{max.intron}{an integer scalar indicating the maximum distance between exons}
\item{key.field}{a character scalar specifying the keytype for name extraction}
\item{name.field}{a character scalar or vector specifying the column(s) to use as the gene name}
\item{ignore.strand}{a logical scalar indicating whether strandedness in \code{incoming} should be ignored}
}

\value{
If \code{incoming} is not provided, a \code{GRanges} object will be returned containing ranges for the exons, promoters and gene bodies. 
Gene keys (e.g., Entrez IDs) are stored as row names.
Gene symbols, exon numbers and internal groupings (for exons of genes with multiple genomic locations) are also stored as metadata.

If \code{incoming} is a \code{GRanges} object, a list will be returned with
\code{overlap}, \code{left} and \code{right} elements. Each element is a
character vector of length equal to the number of ranges in \code{incoming}.
Each non-empty string records the gene symbol, the overlapped exons and the
strand. For \code{left} and \code{right}, the gap between the range and the
annotated feature is also included.
}

\details{
This function adds exon-based annotations to a given set of genomic regions, in the form of compact character strings specifying the features overlapping and flanking each region.
The aim is to determine the genic context of empirically identified regions. 
This allows some basic biological interpretation of binding/marking in those regions. 
All neighboring genes within a specified range are reported, rather than just the closest gene to the region.
If a region in \code{incoming} is stranded and \code{ignore.strand=FALSE}, annotated features will only be reported if they lie on the same strand as that region.

If \code{incoming} is missing, then the annotation will be provided directly to the user in the form of a \code{GRanges} object. 
This may be more useful when further work on the annotation is required. 
Exon numbers are provided in the metadata with promoters and gene bodies labelled as 0 and -1, respectively.
Overlaps to introns can be identified by finding those regions that overlap with gene bodies but not with any of the corresponding exons.
}

\section{Explanation of fields}{
For annotated features overlapping a region, the character string in the \code{overlap} output vector will be of the form \code{GENE|EXONS|STRAND}. 
\code{GENE} is the gene symbol by default, but reverts to \code{<XXX>} if no symbol is defined for a gene with the Entrez ID \code{XXX}.
The \code{EXONS} indicate the exon or range of exons that are overlapped. 
The \code{STRAND} is, obviously, the strand on which the gene is coded. 
For annotated regions flanking the region within a distance of \code{dist}, the character string in the \code{left} or \code{right} output vectors will have an additional \code{[DIST]} value. 
This represents the gap between the edge of the region and the closest exon for that gene.

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
}

\section{Other options}{
The \code{max.intron} value is necessary to deal with genes that have ambiguous locations on the genome. 
If a gene has exons on different chromosomes, its location is uncertain and the gene is partitioned into two sets of exons for separate processing. 
However, this is less obvious when the ambiguous locations belong to the same chromosome. 
The \code{max.intron} value protects against excessively large genes that may occur from considering those locations as a single transcriptional unit.
Exons are partitioned into two (or more) internal groupings for further processing.

The default settings for \code{key.field} and \code{name.field} will work for human and mouse genomes, but may not work for other organisms.
The \code{key.field} should refer to the key type in the \code{OrgDb} object, and also correspond to the \code{GENEID} of the \code{TxDb} object.
For example, in S. cerevisiae, \code{key.field} is set to \code{"ORF"} while \code{name.field} is set to \code{"GENENAME"}.
If multiple entries are supplied in \code{name.field}, the value of \code{GENE} is defined as a semicolon-separated list of each of those entries.
}

\author{Aaron Lun}

\examples{ 
require(org.Mm.eg.db)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)

current <- readRDS(system.file("exdata", "exrange.rds", package="csaw"))
output <- detailRanges(current, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene) 
head(output$overlap)
head(output$right)
head(output$left)

detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db)

\dontrun{
output <- detailRanges(current, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
    orgdb=org.Mm.eg.db, name.field=c("ENTREZID"))
head(output$overlap)

output <- detailRanges(current, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
    orgdb=org.Mm.eg.db, name.field=c("SYMBOL", "ENTREZID"))
head(output$overlap)
}
}

\keyword{annotation}
