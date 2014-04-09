###################################################################################################
# This script tests the PET-analysis capabilities of the 'csaw' package.

suppressPackageStartupMessages(library(csaw))

source("simsam.R")
dir<-"pet-test"
dir.create(dir)

compstat<-function (nreads, singles, chromosomes) {
    # Seeding all reads.
	names<-paste('x', rep(1:nreads, 2), sep=".");
	chrs<-sample(length(chromosomes), length(names), replace=TRUE)
	pos<-rep(0, length(names));

	# Assigning positions to all of them.
	for (i in 1:length(chromosomes)) {
		current<-chrs==i;
		pos[current]<-round(runif(sum(current), 1, chromosomes[i]))
	}

    # Throwing them into the SAM file generator. 
	str<-rbinom(nreads*2, 1, 0.5)==1; 
	outname<-file.path(dir, "out")
	rlen<-10;
	reversi<-c(nreads+1:nreads, 1:nreads)
	out<-simsam(outname, names(chromosomes)[chrs], pos, str, chromosomes, names=names, len=rlen,
			is.paired=TRUE, is.first=c(rep(TRUE, nreads), rep(FALSE, nreads)), 
			mate.chr=names(chromosomes)[chrs][reversi], mate.pos=pos[reversi], mate.str=str[reversi])

	# Adding some singles.
	if (singles) {
		snames<-paste('y', 1:singles, sep=".");
		schrs<-sample(length(chromosomes), singles, replace=TRUE)
		spos<-rep(0, singles);
		for (i in 1:length(chromosomes)) {
			scurrent<-schrs==i;
			spos[scurrent]<-round(runif(sum(scurrent), 1, chromosomes[i]))
		}
		sstr<-rnbinom(singles, 1, 0.5)==1

		tempname<-file.path(dir, "temp")
		out2<-simsam(tempname, names(chromosomes)[schrs], spos, sstr, chromosomes, names=snames, len=rlen)
		more.temp<-file.path(dir, "temp2")
		out<-mergeBam(c(out, out2), more.temp, indexDestination=TRUE, overwrite=TRUE)
		file.rename(more.temp, out)
	}

	stuff<-getPETSizes(out);
	if (stuff$diagnostics[["single"]]!=singles) { 
		stop("mismatch in number of singles");
	} else if (stuff$diagnostics[["total"]]!=nreads*2+singles) {
		stop("mismatch in total number of reads");
	} 

    # Checking whether the sizes make sense.
	pos1<-pos[1:nreads]
	pos2<-pos[nreads+1:nreads]
	chr1<-chrs[1:nreads]
	chr2<-chrs[nreads+1:nreads]
	str1<-str[1:nreads]
	str2<-str[nreads+1:nreads]
	valid<-(chr1==chr2 & str1!=str2 & ifelse(str1, pos1 <= pos2, pos2 <= pos1))

	print(stuff$diagnostics)	
	if (sum(chr1!=chr2)!=stuff$diagnostics[["inter.chr"]]) { stop("mismatch in interchromosomal PETs"); }
    
   	pos1[!str1]<-pmin(pos1[!str1]+rlen, chromosomes[chr1][!str1]+1)
	pos2[!str2]<-pmin(pos2[!str2]+rlen, chromosomes[chr2][!str2]+1)
   	sizes<-abs(pos1-pos2)[valid]
	stopifnot(length(sizes)==length(stuff$sizes));
	if (any(sort(sizes)!=sort(stuff$sizes))) { stop("mismatch in sizes"); }
	
	head(sizes);
}

# Now running through a gamut of tests.

set.seed(1002)

compstat(500, 0, c(chrA=1000, chrB=2000))

compstat(1000, 0, c(chrA=1000, chrB=2000))

compstat(2000, 0, c(chrA=1000, chrB=2000))

compstat(5000, 0, c(chrA=1000, chrB=2000))

# Again, with some singles.

compstat(500, 20, c(chrA=1000, chrB=2000))

compstat(1000, 100, c(chrA=1000, chrB=2000))

compstat(2000, 100, c(chrA=1000, chrB=2000))

compstat(5000, 100, c(chrA=1000, chrB=2000))

# Going pathologically low to force the PET algorithm to process end-of-chromosome catches and stacked F/R reads.
	
compstat(1000, 0, c(chrA=50))

compstat(1000, 0, c(chrA=10))

compstat(2000, 0, c(chrA=10))

compstat(5000, 0, c(chrA=10))

###################################################################################################
# We then write a function to check the extraction procedure.	
	
checkcount<-function (nreads, chromosomes, spacing=50, max.frag=500, left=0, right=0, filter=-1, ext=100) {
	stuff<-file.path(dir, paste("x", 1:2, sep=""));
	all.ranges<-list()
	first.is.left <- list()

	for (x in 1:length(stuff)) {
    	# Seeding all reads.
		names<-paste('x', rep(1:nreads, 2), sep=".");
		chrs<-sample(length(chromosomes), length(names), replace=TRUE)
		pos<-rep(0, length(names));

		# Assigning positions to all of them.
		for (i in 1:length(chromosomes)) {
			current<-chrs==i;
			pos[current]<-round(runif(sum(current), 1, chromosomes[i]))
		}

    	# Throwing them into the SAM file generator. 
		str<-rbinom(length(names), 1, 0.5)==1; 
		rlen<-10;
		reversi<-c(1:nreads+nreads, 1:nreads)
		out<-simsam(stuff[x], names(chromosomes)[chrs], pos, str, chromosomes, names, 
				is.first=c(rep(TRUE, nreads), rep(FALSE, nreads)), is.paired=TRUE, 
				mate.chr=names(chromosomes)[chrs][reversi], mate.pos=pos[reversi], 
				mate.str=str[reversi], len=rlen);

    	# Extracting all valid hits.
		pos1<-pos[1:nreads]
		pos2<-pos[nreads+1:nreads]
		chr1<-chrs[1:nreads]
		chr2<-chrs[nreads+1:nreads]
		str1<-str[1:nreads]
		str2<-str[nreads+1:nreads]
		valid<-(chr1==chr2 & str1!=str2 & ifelse(str1, pos1 <= pos2, pos2 <= pos1))
    	
   		pos1[!str1]<-pmin(pos1[!str1]+rlen, chromosomes[chr1][!str1]+1)
		pos2[!str2]<-pmin(pos2[!str2]+rlen, chromosomes[chr2][!str2]+1)
   		sizes<-abs(pos1-pos2)[valid]
		vpos<-ifelse(str1, pos1, pos2)[valid]
		vchrs<-chr1[valid]
	
		# Filtering by the maximum size.	
		is.spaced<-sizes>max.frag;
		sizes<-sizes[!is.spaced];
		vpos<-vpos[!is.spaced];
	 	vchrs<-vchrs[!is.spaced];

    	# Compiling them into counts for each chromosome.
		all.ranges[[x]]<-GRanges(names(chromosomes)[vchrs], IRanges(vpos, vpos+sizes-1))
		first.is.left[[x]] <- str1[valid][!is.spaced]
	}

    # Collating all counts.
	fnames <- paste0(stuff, ".bam")
	x<-windowCounts(fnames, spacing=spacing, max.frag=max.frag, left=left, right=right, pet="both", filter=filter)
	out<-matrix(0L, length(x$region), length(stuff))
	totals<-integer(length(stuff))
	for (i in 1:length(stuff)) {
        current<-findOverlaps(x$region, all.ranges[[i]])
        out[,i]<-tabulate(queryHits(current), nbins=length(x$region))
		totals[i]<-length(all.ranges[[i]])
	}

	if (!identical(out, x$counts)) { stop('mismatches in counts for paired data') }
	if (!identical(totals, x$totals)) { stop("mismatches in totals for paired data") }

	# Also checking forward and reverse counts.
	for (mode in c("first", "second")) { 
		bravo <- windowCounts(fnames, spacing=spacing, ext=ext, left=left, right=right, pet=mode, filter=filter)
		out<-matrix(0L, length(x$region), length(stuff))
		totals<-integer(length(stuff))
		for (i in 1:length(stuff)) {
			cur.ranges <- all.ranges[[i]]
			if (mode=="first") {
				cur.stat <- first.is.left[[i]]
			} else {
				cur.stat <- !first.is.left[[i]]
			}
			start(cur.ranges) <- ifelse(cur.stat, start(cur.ranges), end(cur.ranges)-ext+1L)
			end(cur.ranges) <- ifelse(cur.stat, start(cur.ranges)+ext-1L, end(cur.ranges))
        	current<-findOverlaps(x$region, all.ranges[[i]])
       		out[,i]<-tabulate(queryHits(current), nbins=length(x$region))
			totals[i]<-length(all.ranges[[i]])
		}
		if (!identical(out, x$counts)) { stop('mismatches in counts for single coercion') }
		if (!identical(totals, x$totals)) { stop("mismatches in totals for single coercion") }
	}
		

	return(x$region)
}

# Running through a bunch of tests.

checkcount(1000, c(chrA=1000, chrB=2000), spacing=20)

checkcount(2000, c(chrA=1000, chrB=2000), spacing=50)

checkcount(5000, c(chrA=1000, chrB=2000), spacing=25)

checkcount(5000, c(chrA=1000, chrB=2000), spacing=25, max.frag=100)

# Checking out restrictions on the max size.

checkcount(1000, c(chrA=1000, chrB=2000), spacing=50, right=0)

checkcount(1000, c(chrA=1000, chrB=2000), spacing=100, right=20)

checkcount(2000, c(chrA=1000, chrB=2000), spacing=100, right=10, max.frag=200)

# Checking out window extension details.

checkcount(1000, c(chrA=1000, chrB=2000), spacing=30, right=100)

checkcount(1000, c(chrA=1000, chrB=2000), spacing=20, left=100)

checkcount(2000, c(chrA=1000, chrB=2000), spacing=15, right=20, left=20)

checkcount(2000, c(chrA=1000, chrB=2000), spacing=15, right=10, left=10, max.frag=200)
	
# Checking out read extension for singles.

checkcount(1000, c(chrA=1000, chrB=2000), spacing=20, ext=100)

checkcount(2000, c(chrA=1000, chrB=2000), spacing=50, ext=50)

checkcount(5000, c(chrA=1000, chrB=2000), spacing=25, ext=20)

checkcount(5000, c(chrA=1000, chrB=2000), spacing=25, ext=200)
	
###################################################################################################
# Cleaning up.

unlink(dir, recursive=TRUE);

###################################################################################################
# End.


