###################################################################################################
# This is a function which just simulates a SAM file and compresses it. It takes a set of names,
# positions and strands, and just fills in the rest.

simsam<-function(f.out, pos.chr, pos.pos, strands, chromosomes, mapq=199, is.dup=NULL, names=NULL, 
		is.first=NULL, is.paired=FALSE, mate.chr=NULL, mate.pos=NULL, mate.str=NULL, len=10) {
	samFile<-paste(f.out, ".sam", sep="")
	out<-file(samFile, open="w");
	for (chr in names(chromosomes)) {
		write(c("@SQ", paste("SN:", chr, sep=""), paste("LN:", chromosomes[[chr]], sep="")), 
				ncolumns=3, file=out, sep="\t");
	}
	cigar<-paste(len, "M", sep="");
	seq<-paste(rep("N", len), collapse="");
	qual<-paste(rep(".", len), collapse="");
	if (is.null(names)) { names<-rep("x", length(pos.pos)); }
	flags<-ifelse(strands, 0, 16)  
	if (length(is.dup)) { flags <- flags + ifelse(is.dup, 1024, 0) }

	# Setting a whole host of paired end information.
	if (is.null(mate.pos)) { mate.pos<-0 }
	if (is.null(mate.chr)) { mate.chr<-"*" }
	isize<-0
	if (is.paired) {
		flags<-flags+1;
        stopifnot(!is.null(is.first)) 
		flags<-flags+ifelse(is.first, 64, 128); 
		flags<-flags+ifelse(mate.str, 0, 32);
		stopifnot(!is.null(mate.pos))
		stopifnot(!is.null(mate.chr))

		mate.chr.size<-chromosomes[match(mate.chr, names(chromosomes))]+1L
		mate.chr<-ifelse(mate.chr==pos.chr, "=", mate.chr)
		isize<-pmin(pmax(mate.pos, pos.pos)+len, mate.chr.size)-pmin(mate.pos, pos.pos)
		isize[mate.pos < pos.pos]<-isize[mate.pos < pos.pos]*-1
		isize[mate.chr!="="]<-0
	}

	mapq <- rep(mapq, length.out=length(names))
	stuff<-data.frame(names, flags, pos.chr, pos.pos, mapq, cigar, mate.chr, mate.pos, isize, seq, qual);
	write.table(file=out, stuff, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
	close(out);

    # We sort and compress to BAM as well.
	suppressPackageStartupMessages(require(Rsamtools))
	tempName<-paste(f.out, "_temp", sep="")
	tempName<-asBam(samFile, destination=tempName, overwrite=TRUE, indexDestination=FALSE);
    newName<-sortBam(tempName, destination=f.out);
	indexBam(newName);
	unlink(tempName);
	return(newName)
}

###################################################################################################

regen <- function(nreads, chromos, outfname) {
	pos.chr<-sample(length(chromos), nreads, replace=TRUE)
	pos.pos<-rep(0, nreads)
	str<-rep(0, nreads)
	for (i in 1:length(chromos)) {
		current<-pos.chr==i
		pos.pos[current]<-round(runif(sum(current), 1, chromos[i]))
		str[current]<-(rbinom(sum(current), 1, 0.5)==1)
	}
	isdup <- rbinom(nreads, 1, 0.8)==0L
    mapq <- round(runif(nreads, 50, 199))
	simsam(outfname, names(chromos)[pos.chr], pos.pos, str, chromos, is.dup=isdup, mapq=mapq)
}

###################################################################################################

makeDiscard <- function(ndisc, sizeof, chromos) {
	chosen <- sample(length(chromos), ndisc, replace=TRUE)
	chosen.pos <- runif(ndisc, 1, chromos[chosen]-sizeof)
	reduce(GRanges(names(chromos)[chosen], IRanges(chosen.pos, chosen.pos+sizeof)))
}

###################################################################################################

generateWindows <- function(chrs, nwin, winsize) {
	allregs<-GRanges()
	for (x in names(chrs)) {
		max.step<-floor(chrs[[x]]/nwin)
		stopifnot(max.step >= 1)
		pos<-cumsum(round(runif(nwin, 1, max.step)))
		suppressWarnings(allregs<-c(allregs, GRanges(x, IRanges(pos, 
			pmin(chrs[[x]], pos+winsize-1L)))))
	}
	total.n<-nwin*length(chrs)
	return(allregs)
}

###################################################################################################

