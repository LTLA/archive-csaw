dumpPE <- function(bam.file, prefix, param=readParam(pe="both"), overwrite=FALSE)
# This function extracts all relevant data, and then dumps it out into another
# SAM file. The idea is to save some time by processing paired-end data once,
# such that we only have to load the position and insert size for each entry.
# We can avoid name loading and matching, which soaks up the most time.
#
# written by Aaron Lun
# created 14 February 2015
# last modified 22 July 2015
{
	if (param$pe!="both") { stop("paired-end inputs required") }
	extracted.chrs <- .activeChrs(bam.file, param$restrict)
	param <- reform(param, fast.pe=FALSE) # So that it properly extracts reads.
	
	# Storing the chromosome lengths.
	ofile <- tempfile(tmpdir='.')
   	on.exit({ if (file.exists(ofile)) { unlink(ofile) } })	
	ohandle <- file(ofile, open='w')
	utils::write.table(file=ohandle, data.frame("@SQ", paste0("SN:", names(extracted.chrs)),
 		paste0("LN:", extracted.chrs)), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)	
	utils::write.table(file=ohandle, data.frame("@PG", "ID:dumped", dump.prog.name),
		sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)	

	# Paired, properly paired, forward strand, mate is reverse strand.
	active.flag <- 0x1 + 0x2 + 0x20 + 0x40 
	counter <- 0L
	all.seq <- all.qual <- character(0)

	for (i in seq_along(extracted.chrs)) {
		chr <- names(extracted.chrs)[i]
		where <- GRanges(chr, IRanges(1L, extracted.chrs[i]))
		out <- .getPairedEnd(bam.file, where=where, param=param, with.reads=TRUE)

		# Setting up width data.
		all.cig <- all.end <- all.names <- list()
		all.widths <- list()
		index <- 1L
		nl <- length(out$left$pos)
		if (nl) { 
			all.cig[[index]] <- paste0(out$left$qwidth, "M")
			all.end[[index]] <- out$right$pos		
			all.names[[index]] <- paste0("paired:", counter+seq_len(nl), ":", out$right$qwidth)
			all.widths[[index]] <- out$left$qwidth
			index <- index + 1L
			counter <- counter + nl
		}
		nr <- length(out$rescued$pos)
		if (nr) { 
			all.cig[[index]] <- paste0(out$rescued$qwidth, "M")
			all.end[[index]] <- out$rescued$pos
			all.names[[index]] <- paste0("rescued:", counter+seq_len(nr), ":", out$rescued$strand)
			all.widths[[index]] <- out$rescued$qwidth
			counter <- counter + nr
		}

		if (length(out$pos)) {
			all.cig <- unlist(all.cig)
			all.end <- unlist(all.end)
			all.names <- unlist(all.names)
			all.widths <- unlist(all.widths)

			# Setting up sequence, quality strings.
			ref.widths <- which(tabulate(all.widths)!=0L)
			if (any(is.na(all.seq[ref.widths]))) {
				for (w in ref.widths) {
					all.seq[w] <- paste(rep("N", w), collapse="")
					all.qual[w] <- paste(rep("#", w), collapse="")
				}
			}	

			# Shorten rescued fragment if it runs off the front (negative positions are not stored).
			subzero <- out$pos <= 0L 
			if (any(subzero)) { 
				ending <- out$pos[subzero] + out$size[subzero] - 1L
				out$pos[subzero] <- 1L
				out$size[subzero] <- ending 
			}

			utils::write.table(file=ohandle, data.frame(all.names, active.flag, chr, out$pos, 255, 
				all.cig, "=", all.end, out$size, all.seq[all.widths], all.qual[all.widths]),
				sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
		}	
	}

	close(ohandle)
	output <- suppressWarnings(asBam(ofile, prefix, indexDestination=TRUE, overwrite=overwrite))
	return(invisible(output))
}

# Assorted helper functions for use elsewhere.

dump.prog.name <- "PN:csaw::dumpPE"
.isDumpedBam <- function(bam.file) {
	return(dump.prog.name %in% scanBamHeader(bam.file, what="text")[[1]]$text[["@PG"]])
}

.isRescuedRead <- function(qnames) {
	return(grepl("^rescued:", qnames))
}

.getMateWidth <- function(qnames) { 
	return(as.integer(sub("^.*:", "", qnames)))
}

.getRescueStr <- function(qnames) {
	return(sub("^.*:", "", qnames))
}

