if [[ $# -eq 0 ]]
then
	RCMD=R
else
	RCMD=$1
fi
shift

trap exit SIGINT SIGTERM
for x in `ls test-*.R`
do
	echo Running $x...
	$RCMD CMD BATCH $@ --no-save --no-restore $x ${x}out
	$RCMD CMD Rdiff ${x}out ${x}out.save
done
