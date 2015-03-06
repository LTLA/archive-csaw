dumpPE <- function(bam.file, prefix, param=readParam(pe="both"), overwrite=FALSE)
# This function extracts all relevant data, and then dumps it out into another
# SAM file. The idea is to save some time by processing paired-end data once,
# such that we only have to load the position and insert size for each entry.
# We can avoid name loading and matching, which soaks up the most time.
#
# written by Aaron Lun
# created 14 February 2015
{
	if (param$pe!="both") { stop("paired-end inputs required") }
	extracted.chrs <- .activeChrs(bam.file, param$restrict)
	param <- reform(param, fast.pe=FALSE) # So that it properly extracts reads.
	
	# Storing the chromosome lengths.
	ofile <- tempfile(tmpdir='.')
   	on.exit({ if (file.exists(ofile)) { unlink(ofile) } })	
	ohandle <- file(ofile, open='w')
	write.table(file=ohandle, data.frame("@SQ", paste0("SN:", names(extracted.chrs)),
 		paste0("LN:", extracted.chrs)), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)	

	counter <- 0L
	active.flag <- 0x1 + 0x2 + 0x20 + 0x40 # Paired, properly paired, forward strand, mate is reverse strand.

    for (i in 1:length(extracted.chrs)) {
		chr <- names(extracted.chrs)[i]
		where <- GRanges(chr, IRanges(1L, extracted.chrs[i]))

        if (.rescueMe(param)) { 
			out <- .rescuePE(bam.file, where=where, param=param)
		} else {
			out <- .extractPE(bam.file, where=where, param=param)
		}

		if (length(out$pos)) {
			fids <- paste0("f", counter + 1:length(out$pos))
			fullsize <- out$size
			ending <- out$pos + out$size - 1L
			subzero <- out$pos <= 0L # Shorten rescued fragment if it runs off the front (negative positions are not stored).
			if (any(subzero)) { 
				out$pos[subzero] <- 1L
				fullsize[subzero] <- ending[subzero] 
			}

			write.table(file=ohandle, data.frame(fids, active.flag, chr, out$pos, 255, "1M", "=", 
				ending, fullsize, "N", "#"), sep="\t", quote=FALSE, col.names=FALSE, 
				row.names=FALSE, append=TRUE)
		}	
		counter <- counter + length(out$pos)
	}

	close(ohandle)
	output <- suppressWarnings(asBam(ofile, prefix, indexDestination=TRUE, overwrite=overwrite))
	return(invisible(output))
}
