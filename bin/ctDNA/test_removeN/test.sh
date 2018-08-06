
# test.sh 20180701 1N >20180701.1N.sh 

time=$1
N=$2

outdir=$time.$N

if [ -d $outdir ] ;then
	`rm -r $outdir`
	`mkdir $outdir`
else
	`mkdir $outdir`
fi
perl_b=/data/bioit/biodata/wugz/small_work/zenghuapingWork/depthCorrect/complex_snp/test1/test1/ctDNA/variant/call_snp.b$N.pl
perl_r=/data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/variant/realign_indel.pl
perl_c=/data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/variant/call_snp.pl
perl_f=/data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/variant/fusion_detect.pl
perl_a=/data/bioit/biodata/zenghp/bin/ctDNA_dev/ctDNA/get_hotspot_detect.pl
perl_a=/data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/get_hotspot_detect.pl
perl_g=/data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/tma_genotyping/cffdna_estimation_and_fetal_genotyping.pl

#samplesheet=/data/bioit/biodata/duyp/Project/Thalassemia/$time\-1.txt
samplesheet=/data/bioit/biodata/duyp/Project/Thalassemia/$3.txt


less $samplesheet |awk '{
print "perl '$perl_b' -i /data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_'$time'/05.variant_detect/"$1".GATK_realign.bam -o '$outdir'/"$1".rN";
print "samtools sort '$outdir'/"$1".rN '$outdir'/"$1".rN";
print "samtools index '$outdir'/"$1".rN.bam";
print "perl '$perl_r' -i '$outdir'/"$1".rN.bam  -o '$outdir'/"$1".rN.new.bam";
print "perl '$perl_c' -i '$outdir'/"$1".rN.new.bam -rp "$6" -o '$outdir'/"$1".vcf";
print "perl '$perl_f' -i '$outdir'/"$1".rN.bam  -k "$1" -od '$outdir'/";
}'

echo "perl $perl_a -id $outdir -is $time.txt -k Total -od $outdir"
echo "perl $perl_g -id $outdir/05.hotspots_result/ -i $time.txt -k Total -od $outdir -e"
