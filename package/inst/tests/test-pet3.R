###################################################################################################
# This script tests the PE-analysis capabilities of the 'csaw' package.

suppressWarnings(suppressPackageStartupMessages(library(csaw)))

source("simsam.R")
dir<-"pe-test"
dir.create(dir)
options(width=120)

checkcount<-function (npairs, nsingles, chromosomes, spacing=50, max.frag=500, left=0, right=0, filter=-1, ext=100) {
	stuff<-file.path(dir, paste("x", 1:2, sep=""))
	firsts <- seconds <- singles <- list()
	rlen <- 10L
	storage.mode(chromosomes) <- "integer"

	for (x in 1:length(stuff)) {
    	# Seeding all reads.
		names<-paste('x', rep(1:npairs, 2), sep=".");
		chrs<-sample(length(chromosomes), length(names), replace=TRUE)
		pos<-rep(0, length(names));

		# Assigning positions to all of them.
		for (i in 1:length(chromosomes)) {
			current<-chrs==i;
			pos[current]<-round(runif(sum(current), 1, chromosomes[i]))
		}

    	# Throwing them into the SAM file generator. 
		str<-rbinom(npairs*2, 1, 0.5)==1L 
		reversi<-c(1:npairs+npairs, 1:npairs)
		mapq <- as.integer(runif(npairs*2, 50, 199))
		dup <- rbinom(npairs*2, 1, 0.8)==0L
		current.chrs <- names(chromosomes)[chrs]
		out<-simsam(stuff[x], current.chrs, pos, str, chromosomes, names=names, 
				is.first=c(rep(TRUE, npairs), rep(FALSE, npairs)), is.paired=TRUE, 
				mate.chr=names(chromosomes)[chrs][reversi], mate.pos=pos[reversi], 
				mate.str=str[reversi], len=rlen, mapq=mapq, is.dup=dup)

		everything <- GRanges(current.chrs, IRanges(pos, pos+rlen-1), strand=!str)
		everything$mapq <- mapq
		everything$dup <- dup
		firsts[[x]] <- everything[1:npairs]
		seconds[[x]] <- everything[1:npairs + npairs]

		# Adding singles.
		if (nsingles) {
			snames <- paste('y', 1:nsingles, sep=".")
            schrs <- sample(length(chromosomes), nsingles, replace=TRUE)
		    spos <- numeric(nsingles)
		    for (i in 1:length(chromosomes)) {
				scurrent<-schrs==i;
				spos[scurrent]<-round(runif(sum(scurrent), 1, chromosomes[i]))
			}
			sstr<-rnbinom(nsingles, 1, 0.5)==1L
				
			mapq <- as.integer(runif(nsingles, 50, 199))
			dup <- rbinom(nsingles, 1, 0.8)==0L
			single.chrs <- names(chromosomes)[schrs]
			singles[[x]] <- GRanges(single.chrs, IRanges(spos, spos+rlen-1L), strand=!sstr)
			singles[[x]]$mapq <- mapq
			singles[[x]]$dup <- dup
			
			tempname<-file.path(dir, "temp")
			out2<-simsam(tempname, single.chrs, spos, sstr, chromosomes, names=snames, len=rlen, mapq=mapq, is.dup=dup)
			more.temp<-file.path(dir, "temp2")
			file.rename(out, more.temp)
			out<-mergeBam(c(more.temp, out2), out, indexDestination=TRUE, overwrite=TRUE)
		}
	}
	
	# Looping through a number of possible filters.
	discard <- GRanges()
	restrict <- NULL
	fnames <- paste0(stuff, ".bam")
	for (filter in 1:4) {
		if (filter==1L) {
			dedup <- TRUE
			minq <- 100L
		} else if (filter==2L) {
			dedup <- FALSE
			minq <- 0L
		} else if (filter==3L) {
			discard <- makeDiscard(10, 50, chromosomes)
		} else {
			discard <- GRanges()
			restrict <- "chrA"
		}

    	# Looping through a number of possible extraction scenarios.
		for (mode in 1:2) {
			if (mode==1L) {
				pe <- "both"
			} else if (mode==2L) {
				pe <- "first"
			} else if (mode==3L) {
                pe <- "second"
            }

			# Loading windowCounts.
			rpam <- readParam(pe=pe, max.frag=max.frag, 
				discard=discard, minq=minq, dedup=dedup, restrict=restrict)
			x <- windowCounts(fnames, spacing=spacing, ext=ext, shift=left, 
				width=right+left+1, filter=0, param=rpam)

			counts <- matrix(0L, nrow=nrow(x), ncol=length(fnames))
			totals <- integer(length(fnames))
			for (lib in 1:length(fnames)) { 
				pos1 <- start(firsts[[lib]])
				chr1 <- as.character(seqnames(firsts[[lib]]))
				str1 <- as.logical(strand(firsts[[lib]])=="+")
				pos2 <- start(seconds[[lib]])
				chr2 <- as.character(seqnames(seconds[[lib]]))
				str2 <- as.logical(strand(seconds[[lib]])=="+")

				valid <- chr1==chr2 & str1!=str2 & ifelse(str1, pos1 <= pos2, pos2 <= pos1)
		   		pos1[!str1] <- pos1[!str1]+rlen
				pos2[!str2] <- pos2[!str2]+rlen
   				sizes<-abs(pos1-pos2)

				# Checking which ones are lost.
				keep1 <- (!dedup | !firsts[[lib]]$dup) & firsts[[lib]]$mapq >= minq
				keep2 <- (!dedup | !seconds[[lib]]$dup) & seconds[[lib]]$mapq >= minq
				if (length(discard)) { 
					keep1 <- keep1 & !overlapsAny(firsts[[lib]], discard, type="within")
					keep2 <- keep2 & !overlapsAny(seconds[[lib]], discard, type="within")
				} 
				if (!is.null(restrict)) { 
					keep1 <- keep1 & chr1 %in% restrict
					keep2 <- keep2 & chr2 %in% restrict
				}
				paired <- keep1 & keep2

				# Checking singles.
				if (nsingles) {
                    schr <- as.character(seqnames(singles[[lib]]))
					skeep <- (!dedup | !singles[[lib]]$dup) & singles[[lib]]$mapq >= minq
					if (!is.null(discard)) { skeep <- skeep & !overlapsAny(singles[[lib]], discard, type="within") } 
					if (!is.null(restrict)) { skeep <- skeep & schr %in% restrict }
				} else {
                    schr <- character(0) 
					skeep <- NULL 
				}
	
				###############################################
				# Checking diagnostics.
				if (mode==1L) { 
		        	stuff<-getPESizes(fnames[lib], readParam(pe="both", minq=minq, dedup=dedup, restrict=restrict, discard=discard))
					if (stuff$diagnostics[["single"]]!=sum(skeep)) { 
						stop("mismatch in number of singles")
					} else if (stuff$diagnostics[["mapped.reads"]]!=sum(keep1)+sum(keep2)+sum(skeep)) {
						stop("mismatch in number of mapped reads")
					}
                    if (stuff$diagnostics[["total.reads"]]!=npairs*2L+nsingles) {
                        stop("mismatch in total number of reads")
                    }
			        if (sum(paired & chr1!=chr2)!=stuff$diagnostics[["inter.chr"]]) { stop("mismatch in interchromosomal PEs") }
			        if (sum(paired & chr1==chr2 & !valid)!=stuff$diagnostics[["unoriented"]]) { stop("mismatch in invalid numbers") }
			        if (sum(keep1!=keep2)!=stuff$diagnostics[["mate.unmapped"]]) { stop("mismatch in unmapped numbers") }
			        if (!identical(sort(sizes[valid&paired]), sort(stuff$sizes))) { stop("mismatch in sizes"); }
					if (lib==1L) { print(stuff$diagnostics) }
				}

				###############################################
				# Now, counting; going through and seeing up the valid paired ones.

				leftpos <- pmin(pos1, pos2)
				valid <- valid & sizes <= max.frag
				if (pe=="both") {
					pairedness <- GRanges(chr1, IRanges(leftpos, leftpos+sizes-1))[valid & paired]
				} else {
					pairedness <- resize(firsts[[lib]][keep1], width=ext)
				}
				counts[,lib] <- countOverlaps(rowRanges(x), pairedness)
				totals[lib] <- length(pairedness)
			}
#			print(c(totals, x$totals))
#			print(which(counts!=x$counts))
#			print(is.integer(counts))
#			print(is.integer(x$counts))
#			print(head(counts))
#			print(head(x$counts))

			curcounts <- assay(x)
			attributes(curcounts)$dimnames <- NULL
			if (!identical(counts, curcounts)) { stop('mismatches in counts for paired data') }
			if (!identical(totals, x$totals)) { stop("mismatches in totals for paired data") }

			# Comparing windowCounts to regionCounts.
			x2 <- regionCounts(fnames, regions=rowRanges(x), ext=ext, param=rpam)
			stopifnot(identical(assay(x), assay(x2)))
			stopifnot(identical(colData(x), colData(x2)))
			stopifnot(identical(rowRanges(x), rowRanges(x2)))

			# Comparing regionCounts to extractReads, using the middle region.
			chosen <- round(nrow(x)/2)
			my.reg <- rowRanges(x)[chosen]
			if (rpam$pe=="both") {
				for (f in 1:length(fnames)) {
					collected <- extractReads(fnames[f], my.reg, param=rpam)
					if (!identical(assay(x)[chosen,f], length(collected))) {
						stop("mismatch in the number of fragments from extractReads")
					}
				}
			} else {
				for (f in 1:length(fnames)) {
					collected <- extractReads(fnames[f], my.reg, param=rpam, ext=ext)
					if (!identical(assay(x)[chosen,f], length(collected))) { 
						stop("mismatch in the number of single reads from extractReads")
					}
				}
			}
		}
	}
	return(rowRanges(x))
}

# Running through a bunch of tests.

set.seed(3429201)
checkcount(1000, 50, c(chrA=1000, chrB=2000), spacing=20)

checkcount(2000, 0, c(chrA=1000, chrB=2000), spacing=50)

checkcount(5000, 25, c(chrA=1000, chrB=2000), spacing=25)

checkcount(5000, 10, c(chrA=1000, chrB=2000), spacing=25, max.frag=100)

# Checking out restrictions on the max size.

checkcount(1000, 10, c(chrA=1000, chrB=2000), spacing=50, right=0)

checkcount(1000, 20, c(chrA=1000, chrB=2000), spacing=100, right=20)

checkcount(2000, 50, c(chrA=1000, chrB=2000), spacing=100, right=10, max.frag=200)

# Checking out window extension details.

checkcount(1000, 100, c(chrA=1000, chrB=2000), spacing=30, right=100)

checkcount(1000, 0, c(chrA=1000, chrB=2000), spacing=20, left=10)

checkcount(2000, 50, c(chrA=1000, chrB=2000), spacing=15, right=-5, left=10)

checkcount(2000, 25, c(chrA=1000, chrB=2000), spacing=15, right=10, left=10, max.frag=200)
	
# Checking out read extension for singles.

checkcount(1000, 0, c(chrA=1000, chrB=2000), spacing=20, ext=100)

checkcount(2000, 50, c(chrA=1000, chrB=2000), spacing=50, ext=50)

checkcount(5000, 10, c(chrA=1000, chrB=2000), spacing=25, ext=20)

checkcount(5000, 20, c(chrA=1000, chrB=2000), spacing=25, ext=200)
	
###################################################################################################
# Checking out behaviour with non-trivial CIGAR strings.

suppressPackageStartupMessages(require(GenomicAlignments))
getFragSizes <- function(positions, cigars, include.clip=TRUE) {
    left.cig <- cigars[1]
    right.cig <- cigars[2]
    left.pos <- positions[1]
    right.pos <- positions[2]

    # Sanity check.
    remaining <- right.pos - left.pos - cigarWidthAlongReferenceSpace(left.cig) + cigarWidthAlongReferenceSpace(right.cig)
    stopifnot(remaining >= 0L)

    left.cig <- cigarToRleList(left.cig)[[1]]
    new.left.pos <- left.pos
    if (include.clip) { 
        new.left.pos <- new.left.pos - ifelse(runValue(left.cig)[1]=="S", runLength(left.cig)[1], 0L)
    }
    new.right.pos <- right.pos + cigarWidthAlongReferenceSpace(right.cig)
    right.cig <- rev(cigarToRleList(right.cig)[[1]])
    if (include.clip) { 
        new.right.pos <- new.right.pos +  ifelse(runValue(right.cig)[1]=="S", runLength(right.cig)[1], 0L)
    }

    return(new.right.pos - new.left.pos)
}

chromosomes <- c(chrA=100, chrB=200)
chr <- c("chrA", "chrA")
positions <- c(5, 10)

for (positions in list(
                       c(5L, 10L),
                       c(6L, 100L),
                       c(10L, 20L)
                       )) {
    for (cigars in list(
                        c("5S5M", "5M5S"),
                        c("2S3M5S", "1S8M1S"),
                        c("7S3M", "10M"),
                        c("10M", "8M2S"),
                        c("5M5S", "2S8M")
                        )) {
        out <- simsam(file.path(dir, "test"), chr, positions, c(TRUE, FALSE), chromosomes, is.first=c(TRUE, FALSE),
                      names=c("x.1", "x.1"), is.paired=TRUE, mate.chr=rev(chr), mate.pos=rev(positions), mate.str=c(FALSE, TRUE), cigar=cigars)
        stopifnot(identical(getPESizes(out)$sizes, getFragSizes(positions, cigars)))
        out2 <- csaw:::.getPairedEnd(out, GRanges("chrA", IRanges(1, chromosomes[1])), param=readParam(pe="both"))
        stopifnot(identical(out2$pos, positions[1]))
        stopifnot(identical(out2$size, getFragSizes(positions, cigars, include.clip=FALSE)))
        cat(sprintf("%i (%s), %i (%s), %i", positions[1], cigars[1], positions[2], cigars[2], getPESizes(out)$sizes), "\n")
    }
}

###################################################################################################
# Cleaning up.

unlink(dir, recursive=TRUE)

###################################################################################################
# End.


