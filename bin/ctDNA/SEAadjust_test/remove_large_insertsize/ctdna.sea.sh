perlfusion=/data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/variant/fusion_detect.pl

time=20180721
outdir=/data/bioit/biodata/wugz/wugz_testout/ctDNA/SEAadjust_test/remove_large_insertsize
timedir=$outdir/$time
if [ -d $timedir ] ;then
	`rm -r $timedir`
	`mkdir $timedir`
else
	`mkdir $timedir`
fi

less $outdir/SEA.$time.txt |awk '{
	print "samtools view -H /data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_'$time'/05.variant_detect/"$1".GATK_realign.bam >'$timedir'/"$1".sam";
	print "samtools view /data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_'$time'/05.variant_detect/"$1".GATK_realign.bam |awk  '"'"'{if($9<=350){print $0}}'"'"' >>'$timedir'/"$1".sam";
	print "samtools view -b -S '$timedir'/"$1".sam >'$timedir'/"$1".bam";
	print "samtools sort '$timedir'/"$1".bam '$timedir'/"$1;
	print "samtools index '$timedir'/"$1".bam";
	print "perl '$perlfusion' -i '$timedir'/"$1".bam  -k "$1" -od '$timedir'/";
	print "rm '$timedir'/"$1".sam '$timedir'/"$1".bam";
}'


