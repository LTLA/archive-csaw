getPESizes <- function(bam.file, param=readParam(pe="both")) 
# This function takes a BAM file and reads it to parse the size of the PE fragments. It then
# returns a vector of sizes which can be plotted for diagnostics. The length of the vector
# will also tell you how many read pairs were considered valid. The total number of reads, the
# number of singletons and the number of interchromosomal pairs is also reported.
# 
# written by Aaron Lun
# a long long time ago
# last modified 15 December 2015
{
	if (param$pe!="both") { stop("paired-end inputs required") }
	extracted.chrs <- .activeChrs(bam.file, param$restrict)
    norm.list <- list()
    totals <- singles <- one.unmapped <- mapped <- unoriented <- 0L
    loose.names.1 <- loose.names.2 <- list()

	for (i in seq_along(extracted.chrs)) {
        cur.chr <- names(extracted.chrs)[i]
        output <- .extractPE(bam.file, GRanges(cur.chr, IRanges(1L, extracted.chrs[i])), param=param, diagnostics=TRUE)
        totals <- totals + output$total
        singles <- singles + output$single

        # Filtering out based on discard.
        relevant <- seqnames(param$discard)==cur.chr
        discard <- param$discard[relevant]

        # For valid read pairs; either go to 'mapped' (and then norm.list), 'one.unmapped', or implicitly unmapped.
        dfkeep <- .discardReads(cur.chr, output$forward[[1]], output$forward[[2]], discard)
        drkeep <- .discardReads(cur.chr, output$reverse[[1]], output$reverse[[2]], discard)
        dkeep <- dfkeep & drkeep 
        mapped <- mapped + sum(dkeep)
        one.unmapped <- one.unmapped + sum(dfkeep!=drkeep)
        all.sizes <- .getFragmentSizes(output$forward, output$reverse)
        norm.list[[i]] <- all.sizes$full[dkeep]

        # For unoriented read pairs; either go to 'mapped' (and the 'unoriented'), 'one.unmapped', or implicitly unmapped.
        ufkeep <- .discardReads(cur.chr, output$ufirst[[1]], output$forward[[2]], discard)
        uskeep <- .discardReads(cur.chr, output$usecond[[1]], output$reverse[[2]], discard)
        ukeep <- ufkeep & uskeep 
        mapped <- mapped + sum(ukeep)
        unoriented <- unoriented + sum(ukeep)
        one.unmapped <- one.unmapped + sum(ufkeep!=uskeep)

        # For lone mappers; either go to 'one.unmapped' or implicitly unmapped.
        omkeep <- .discardReads(cur.chr, output$one.mapped[[1]], output$one.mapped[[2]], discard)
        one.unmapped <- one.unmapped + sum(omkeep)

        # For inter-chromosomals; either store names, or don't.
        ikeep1 <- .discardReads(cur.chr, output$ifirst[[1]], output$ifirst[[2]], discard)    
        loose.names.1[[i]] <- output$ifirst[[3]][ikeep1]
        ikeep2 <- .discardReads(cur.chr, output$isecond[[1]], output$isecond[[2]], discard)    
        loose.names.2[[i]] <- output$isecond[[3]][ikeep2]
    }

	# Checking whether a read is positively matched to a mapped counterpart on another chromosome.
	# If not, then it's just a read in an unmapped pair.
	loose.names.1 <- unlist(loose.names.1)
	loose.names.2 <- unlist(loose.names.2)
	inter.chr <- sum(loose.names.1 %in% loose.names.2)
	one.unmapped <- one.unmapped + length(loose.names.2) + length(loose.names.1) - inter.chr*2L

   	# Returning sizes and some diagnostic data.
	return(list(sizes=unlist(norm.list), diagnostics=c(total.reads=totals, mapped.reads=mapped, 
		single=singles, mate.unmapped=one.unmapped, unoriented=unoriented, inter.chr=inter.chr)))
}

.getFragmentSizes <- function(left, right) {
    clipped.sizes <- right[[1]] - left[[1]] + right[[2]]
    full.sizes <- clipped.sizes + left[[3]] + right[[3]]
    list(clipped=clipped.sizes, full=full.sizes)
}
