###################################################################################################
# This is a function which just simulates a SAM file and compresses it. It takes a set of names,
# positions and strands, and just fills in the rest.

simsam<-function(f.out, pos.chr, pos.pos, strands, chromosomes, names=NULL, is.first=NULL, 
		is.paired=FALSE, mate.chr=NULL, mate.pos=NULL, mate.str=NULL, 
		len=10) {
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

	stuff<-data.frame(names, flags, pos.chr, pos.pos, 199, cigar, mate.chr, mate.pos, isize, seq, qual);
	write.table(file=out, stuff, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
	close(out);

    # We sort and compress to BAM as well.
	require(Rsamtools)
	tempName<-paste(f.out, "_temp", sep="")
	tempName<-asBam(samFile, destination=tempName, overwrite=TRUE, indexDestination=FALSE);
    newName<-sortBam(tempName, destination=f.out);
	indexBam(newName);
	unlink(tempName);
	return(newName)
}

###################################################################################################

