sam=$1
bam=`echo $1 | sed "s/sam$/bam/"`
sorted=`echo $bam | sed "s/.bam$//"`"_sort"
samtools view -bS -o $bam $sam 
samtools sort $bam $sorted
samtools index $sorted".bam"
rm $bam

