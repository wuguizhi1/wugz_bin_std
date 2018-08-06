
time=$1
dirout=/data/bioit/biodata/wugz/wugz_testout/ctDNA/SEAadjust_test/remove_large_insertsize/

less /data/bioit/biodata/duyp/Project/Thalassemia/$time.txt |cut -f 1 |sort -u >$dirout/SEA.$time.txt

#timedir="$dirout/$time/adjust"
#if [ -d $timedir ] ;then
#	`rm -r $timedir`
#	`mkdir $timedir`
#else
#	`mkdir $timedir`
#fi

hotspotdir=/data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$time.2N/05.hotspots_result/
hotspotdir=/data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_$time/statistic/05.hotspots_result/
less $dirout/SEA.$time.txt |awk '{
	print "perl SEA_intersize_pipe.pl -t '$time' -s "$1" -od '$dirout$time'";
	print "less '$hotspotdir'/"$1".hotspots.result |grep -v SEA >'$dirout$time'/"$1".hotspots.result";
	print "less '$dirout$time'/"$1".fusion.final.txt |perl grep_SEA.pl '$dirout$time' "$1" >>'$dirout$time'/"$1".hotspots.result";
	print "Rscript /data/bioit/biodata/wugz/wugz_bin/bin/ctDNA/cffDNAadjust/bin/refine_fetal_genotyping_EM.old.r -i '$dirout$time'/"$1".hotspots.result -o '$dirout$time'/"$1".genotype.txt -d 100 -s 3 -a 2";
}'
	echo "perl InsertSize.pl -i $dirout/SEA.$time.txt -id $dirout$time -c /data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_$time/06.genotyping/cffdna.txt -od $timedir -m 10000 -k SEA"
less /data/bioit/biodata/wugz/wugz_testout/ctDNA/SEAadjust_test/remove_large_insertsize/$time/adjust/SEA.size.txt.adjust.txt |awk '{if($2 ~ /SEA/){ print "Rscript /data/bioit/biodata/wugz/wugz_bin/bin/ctDNA_together/cffDNA_to_geno.onbyone.r -i /data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/'$time'.2N/no_removeN_05.hotspots_result/"$1".hotspots.result -c "$7" -o /data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/'$time'.2N/no_removeN_05.hotspots_result/"$1".genotype.txt" }}'"
