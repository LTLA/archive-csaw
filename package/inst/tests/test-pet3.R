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
		for (mode in 1:3) {
			if (mode==1L) {
				rescue.ext <- NA
				pe <- "both"
			} else if (mode==2L) {
				rescue.ext <- ext
			} else if (mode==3L) {
				pe <- "first"
			}

			# Loading windowCounts.
			rpam <- readParam(pe=pe, rescue.ext=rescue.ext, max.frag=max.frag, 
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
					if (is.null(restrict)) { 
                        if (stuff$diagnostics[["total.reads"]]!=npairs*2L+nsingles) {
                            stop("mismatch in total number of reads")
                        }
                    } else {
                        if (stuff$diagnostics[["total.reads"]]!=sum(chr1%in%restrict) + sum(chr2%in% restrict) + sum(schr%in% restrict)) {
                            stop("mismatch in total number of reads upon restriction")
                        }
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
					if (!is.na(rescue.ext)) {
						# We pick the first if the second is inactive, if paired but interchromosomal, or if paired and intrachromosomal
						# and otherwise invalid and has higher mapping quality.
						better <- firsts[[lib]]$mapq > seconds[[lib]]$mapq
						fcopy <- resize(firsts[[lib]][keep1 & (!keep2 | chr1!=chr2 | (!valid & better))], width=ext)
						fcopy$mapq <- fcopy$dup <- NULL
						scopy <- resize(seconds[[lib]][keep2 & (!keep1 | chr1!=chr2 | (!valid & !better))], width=ext)
						scopy$mapq <- scopy$dup <- NULL
						pairedness <- c(pairedness, fcopy, scopy)
					}
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
# Cleaning up.

unlink(dir, recursive=TRUE)

###################################################################################################
# End.


