getPESizes <- function(bam.file, param=readParam(pe="both")) 
# This function takes a BAM file and reads it to parse the size of the PE fragments. It then
# returns a vector of sizes which can be plotted for diagnostics. The length of the vector
# will also tell you how many read pairs were considered valid. The total number of reads, the
# number of singletons and the number of interchromosomal pairs is also reported.
# 
# written by Aaron Lun
# a long long time ago
# last modified 22 July 2015
{
	if (param$pe!="both") { stop("paired-end inputs required") }
	extracted.chrs <- .activeChrs(bam.file, param$restrict)
    bam.file <- path.expand(bam.file)
    bam.index <- paste0(bam.file, ".bai")
    diagnostics <- integer(5)
    norm.list <- list()
    loose.names.1 <- loose.names.2 <- list()

	for (i in seq_along(extracted.chrs)) {
        out <- .Call(cxx_extract_pair_data, bam.file, bam.index, names(extracted.chrs)[i], param$minq, param$dedup, TRUE)
        if (is.character(out)) { stop(out) }
        diagnostics <- diagnostics + out[[3]]
        loose.names.1[[i]] <- out[[4]][[1]]
        loose.names.2[[i]] <- out[[4]][[2]]       
        norm.list[[i]] <- out[[2]][[1]] - out[[1]][[1]] + out[[2]][[2]]   
    }

	# Checking whether a read is positively matched to a mapped counterpart on another chromosome.
	# If not, then it's just a read in an unmapped pair.
	loose.names.1 <- unlist(loose.names.1)
	loose.names.2 <- unlist(loose.names.2)
	inter.chr <- sum(loose.names.1 %in% loose.names.2)
	extra.unmapped <- length(loose.names.2) + length(loose.names.1) - inter.chr*2L

   	# Returning sizes and some diagnostic data.
	return(list(sizes=unlist(norm.list), diagnostics=c(total.reads=diagnostics[1], mapped.reads=diagnostics[2],
		single=diagnostics[3], mate.unmapped=diagnostics[4]+extra.unmapped, unoriented=diagnostics[5], inter.chr=inter.chr)))
}

##################################

.yieldInterestingBits <- function(reads, clen, max.frag=Inf, diag=FALSE)
# This figures out what the interesting reads are, given a list of 
# read names, flags, positions and qwidths. In particular, reads have
# to be properly paired in an inward conformation, without one read
# overrunning the other (if the read lengths are variable).
#
# written by Aaron Lun
# created 15 April 2014
# last modified 13 May 2015
{ 
	if (max.frag <= 0L) {
		return(list(pos=integer(0), size=integer(0), is.ok=logical(length(reads$flag))))
	}

	# Figuring out the strandedness and the read number.
 	is.forward <- bitwAnd(reads$flag, 0x10) == 0L 
 	is.mate.reverse <- bitwAnd(reads$flag, 0x20) != 0L
	should.be.left <- is.forward & is.mate.reverse
	should.be.right <- !is.forward & !is.mate.reverse
	is.first <- bitwAnd(reads$flag, 0x40) != 0L
	is.second <- bitwAnd(reads$flag, 0x80) != 0L
	stopifnot(all(is.first!=is.second))

	# Matching the reads in each pair so only valid PEs are formed.
	set.left.A <- should.be.left & is.first
	set.right.A <- should.be.right & is.second
	set.left.B <- should.be.left & is.second
	set.right.B <- should.be.right & is.first
	corresponding.1 <- match(reads$qname[set.left.A], reads$qname[set.right.A])
	corresponding.2 <- match(reads$qname[set.left.B], reads$qname[set.right.B])

	hasmatch <- logical(length(reads$flag))
	hasmatch[set.left.A] <- !is.na(corresponding.1)
	hasmatch[set.left.B] <- !is.na(corresponding.2)
	corresponding <- integer(length(reads$flag))
	corresponding[set.left.A] <- which(set.right.A)[corresponding.1]
	corresponding[set.left.B] <- which(set.right.B)[corresponding.2]

	# Selecting the read pairs.
	fpos <- reads$pos[hasmatch]
	fwidth <- reads$qwidth[hasmatch]
	rpos <- reads$pos[corresponding[hasmatch]]
	rwidth <- reads$qwidth[corresponding[hasmatch]]

	# Allowing only valid PEs.
	fend <- pmin(fpos+fwidth, clen+1L)
	rend <- pmin(rpos+rwidth, clen+1L)
	total.size <- rend - fpos
	valid <- fpos <= rpos & fend <= rend & total.size <= max.frag
	fpos <- fpos[valid]
	total.size <- total.size[valid]

	output <- list(pos=fpos, size=total.size)
	if (diag) {
		# Don't bother returning matches of invalids; you'd have to re'match
		# to account for improperly oriented pairs anyway, in getPESizes.
		output$left <- which(hasmatch)[valid]
		output$right <- corresponding[hasmatch][valid]
		is.ok <- logical(length(reads$flag))
		is.ok[output$left] <- TRUE
		is.ok[output$right] <- TRUE
		output$is.ok <- is.ok
	}
	return(output)
}

##################################
