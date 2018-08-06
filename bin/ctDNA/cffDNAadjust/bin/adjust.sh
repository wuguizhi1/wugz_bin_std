time=$1
dir=/data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/cffDNAadjust/

if [ -d $dir/$time ] ;then
	`rm -r $dir/$time`
	`mkdir $dir/$time`
	`mkdir $dir/$time/05.hotspots_result`
	`ln -s /data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$time.2N/cffdna.txt $dir/$time/cffdna.old.txt`
	`ln -s /data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$time.2N/05.hotspots_result/* $dir/$time/05.hotspots_result`
else
	`mkdir $dir/$time`
	`mkdir $dir/$time/05.hotspots_result`
	`ln -s /data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$time.2N/cffdna.txt $dir/$time/cffdna.old.txt`
	`ln -s /data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$time.2N/05.hotspots_result/* $dir/$time/05.hotspots_result`
fi

cmd="perl $dir/bin/cffdna_estimation_and_fetal_genotyping.pl -id $dir/$time/05.hotspots_result/ -i /data/bioit/biodata/duyp/Project/Thalassemia/$time.txt -k Total -od $dir/$time/ -e"
echo $cmd


outdir="$dir/$time/05.hotspots_result"
less old.txt |awk 'print " Rscript /data/bioit/biodata/wugz/wugz_bin_std/bin/ctDNA/cffDNAadjust/bin/refine_fetal_genotyping_EM.r -i '$outdir'"$1".hotspots.result -o '$outdir'"$1".genotype.txt -d 100 -s 3 -a 2 \n perl /data/bioit/biodata/wugz/wugz_bin/bin/ctDNA_together/trans_TF.one.pl -i '$outdir'"$1".genotype.txt -g "$2":"$3":"$4" -o '$outdir'"$1".trueorfalse"' >new.txt

lnsdir="/data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/$time.2N/05.hotspots_result/"
less old.txt |awk 'print " Rscript /data/bioit/biodata/wugz/wugz_bin_std/bin/ctDNA/cffDNAadjust/bin/refine_fetal_genotyping_EM.r -i '$lnsdir'"$1".hotspots.result -o '$lnsdir'"$1".genotype.txt -d 100 -s 3 -a 2 \n perl /data/bioit/biodata/wugz/wugz_bin/bin/ctDNA_together/trans_TF.one.pl -i '$lnsdir'"$1".genotype.txt -g "$2":"$3":"$4" -o '$lnsdir'"$1".trueorfalse"' >>new.txt


