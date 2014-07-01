current <- commandArgs(trailingOnly=TRUE)
rcmd <- "R"
all.opts <- ""
if (length(current)) {
	rcmd <- current[1]
	all.opts <- paste(current[-1], collapse=" ")
}

for (x in list.files(".", pattern="^test.*\\.R$")) {
	cat("Running", x, "...")
	system(sprintf("%s CMD BATCH %s --no-save --no-restore %s %sout", rcmd, all.opts, x, x))
	is.diff <- system(sprintf("%s CMD Rdiff %sout %sout.save", rcmd, x, x), intern=TRUE)
	if (length(is.diff)) { 
		cat('\n')
		cat(is.diff, sep='\n')
		break
	}
	cat(' OK\n')
}
