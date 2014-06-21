for x in `ls test-*.R`
do
	echo Running $x...
	R CMD BATCH --no-save --no-restore $x ${x}out
	R CMD Rdiff ${x}out ${x}out.save
done
