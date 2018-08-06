
removeN=$1
#removeN=20180721_25_merge

remove2Ndir=/data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$removeN.2N/05.hotspots_result/
SEAdir=/data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$removeN.2N/no_removeN_05.hotspots_result/
samplesheet=/data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$removeN.txt
rm $SEAdir/*hotspots*

less /data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_$removeN/statistic/Total.genotype.result |grep SEA |cut -f 1 >$removeN.txt
perl InsertSize.pl -i $removeN.txt -id /data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_$removeN -c /data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_$removeN/06.genotyping/cffdna.txt -od $SEAdir -m 10000 -k SEA >$removeN.o 2>$removeN.e

cd $SEAdir
for i in $(ls *hotspots.result)
do
	#echo "#$remove2Ndir/$i\t#$SEAdir/$i"
	less $remove2Ndir/$i |sed -n '2,$p' |grep -v "SEA" >>$SEAdir/$i
done

cd $remove2Ndir
for i in $(ls *hotspots.result)
do
	if [ ! -f $SEAdir/$i ];then
		`ln -s $remove2Ndir/$i $SEAdir/$i`
	fi
done

#less $SEAdir/SEA.size.txt.adjust.txt |awk '{if($2 ~ /SEA/){ print "Rscript /data/bioit/biodata/wugz/wugz_bin/bin/ctDNA_together/cffDNA_to_geno.onbyone.r -i '$SEAdir'/"$1".hotspots.result -c "$7" -o '$SEAdir'/"$1".genotype.txt" }}'

# only SEA sample
#samplesheetgrep=`less $SEAdir/SEA.size.txt.adjust.txt |awk '{if($2 ~ /SEA/){ print $1 }}' |xargs |sed 's/ /\\\|/g'`
#echo "less $samplesheet |grep \"$samplesheetgrep\" >$SEAdir/samplesheet.txt"
#echo "perl /data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/tma_genotyping/cffdna_estimation_and_fetal_genotyping.pl -id $SEAdir -i $SEAdir/samplesheet.txt -k Total -od $SEAdir -e "

echo "perl /data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/tma_genotyping/cffdna_estimation_and_fetal_genotyping.pl -id $SEAdir -i $samplesheet -k Total -od $SEAdir -e "
perl /data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/tma_genotyping/cffdna_estimation_and_fetal_genotyping.pl -id $SEAdir -i $samplesheet -k Total -od $SEAdir -e 

