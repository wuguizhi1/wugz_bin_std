

inbam=/data/bioit/biodata/wugz/small_work/zenghuapingWork/depthCorrect/company_snp/bin/test.sorted.bam
# perl="perl ../call_snp.pl -i $inbam -o $inbam.vcf"     # 364s

perl="perl call_snp_complex.pl -i $inbam -o $inbam.out 1>$inbam.o 2>$inbam.e &"
echo $perl
#$perl

