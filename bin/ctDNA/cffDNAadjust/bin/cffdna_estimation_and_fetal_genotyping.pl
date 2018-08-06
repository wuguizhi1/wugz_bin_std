#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
require "$Bin/../mylib.pl";
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ==============================================================
# Get Options
# ==============================================================
my $fsample;
my $indir;
my $dOut = "./"; 
my $min_depth = 100;
my $infer_genptype_by_cffDNA;
my $fkey;
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fsample,
	"id:s"=>\$indir,
	"od:s"=>\$dOut,
	"k:s"=>\$fkey,
	"d:i"=>\$min_depth,
	"e"=>\$infer_genptype_by_cffDNA,
) or &USAGE;
&USAGE unless ($fsample and $indir and $fkey);

mkdir $dOut if (defined $dOut && not -d $dOut);
$dOut = Cwd::abs_path($dOut);

#===============================================================
# Main
#===============================================================

#样本类型	备注	数据量(M)	样本编号	基因型
#地贫白细胞	引物panel A10，8173引物测试，CD26样本间测试	5	GX077-1	SEA/N,CD26/N
#地贫白细胞	引物panel A10-T8，8173引物测试，CD26样本间测试	5	GX077-1	SEA/N,CD26/N
#样本类型	备注	数据量(M)	样本编号	基因型	比例
#地贫标准品	TMA-(AB-10S-A10)引物，三种预文库建库筛选，标准法	5	GX070-1/GX019	SEA（-/-）	3%
#地贫标准品	TMA-(AB-10S-A10)引物，三种预文库建库筛选，一管法	5	GX070-1/GX019	SEA（-/-）	3%
#样本类型	备注	数据量(M)	样本编号	胎儿基因型	母亲基因型	父亲基因型	plasma vol (ml)	DNA con.(ng/ul)	prelib  con.(ng/ul)
#地贫检测孕妇血浆	TMA-(AB-10S-A10-2)引物+SNP，真实样本检测，盲测	12	FS292				0.8	0.353	5.9
#地贫检测孕妇血浆	TMA-(AB-10S-A10-2)引物+SNP，真实样本重新检测，低深度	12	NF386	IVS-II-654/N,-28/N	 -28/N	IVS-II-654/N	2	0.452	9.16
my %config;
&read_data_sheet($fsample, \%config); ## @{$adata->{$sample}->{"sample_info"}}=@sample_info;

# write shell for fetal genotyping and fetal fraction estimation using EM algorithm
my $sh = "$dOut/cffdna_genotype.sh";
open (SH,">$sh") or die $!;
foreach my $s (sort {$a cmp $b} keys %config) {
	my @saminfo = @{$config{$s}{"sample_info"}};
	my $fin = "$indir/$s.hotspots.result";
	my $fout = "$dOut/$s.genotype.txt";
	my $cmd;
	next if($config{$s}{"sample_info"}->[0] ne "地贫孕妇血浆" && $config{$s}{"sample_info"}->[0] ne "地贫标准品");
	if($saminfo[0] eq "地贫孕妇血浆"){
		$cmd = "Rscript $Bin/refine_fetal_genotyping_EM.r -i $fin -o $fout -d $min_depth -s 3 -a 2 ";
		$cmd .= "-e " if (defined $infer_genptype_by_cffDNA);
	}elsif($saminfo[0] eq "地贫标准品"){
		my $mix_ratio = ($saminfo[5]*100+0.1483)/0.6597/100;
		$cmd = "Rscript $Bin/cffDNA_to_genotype.r -i $fin -c $mix_ratio -o $fout";
	}
	if($saminfo[4] ne ""){
		my $genotype;
		if($saminfo[0] eq "地贫孕妇血浆"){
			$genotype = join(":", @saminfo[4..6]);
		}elsif($saminfo[0] eq "地贫标准品"){
			my ($spot, @geno)=$saminfo[4]=~/(\S+)\(([-N]+)\/([-N]+)\)/;
			for(my $i=0; $i<@geno; $i++){
				if($geno[$i] eq "-"){
					$geno[$i] = $spot."/N";
				}elsif($geno[$i] eq "--"){
					$geno[$i] = $spot."/".$spot;
				}elsif($geno[$i] eq "N"){
					$geno[$i] = "N/N";
				}else{
					die "Wrong info $saminfo[4] of $s\n";
				}
			}
			$genotype = join(":",$geno[1],$geno[0],"");
		}

		$cmd .= " && perl $Bin/TF_evaluation.pl -i $fout -g $genotype -o $fout";
	}
	print SH $cmd, "\n"; 
}
close (SH) ;

my $r = system("parallel -j 10 < $sh ");
die "infer fetal genotypes failed\n" if ($r);


## merge result
open(O,">$dOut/$fkey.genotype.result") or die $!;
print O "Sample\tMutID\tMutTotal\tDepTotal\tFreqTotal\tfetal.EM\tMaternalGenotype\tFetalGenotype\tRefined.genotype.EM\tp\tAverDepth\tAverDepthHotspot\tCheck\tFetalGiving\tMotherGiving\tFatherGiving\tTrueOrFalse\n";
foreach my $s (sort {$a cmp $b} keys %config) {
	next if($config{$s}{"sample_info"}->[0] ne "地贫孕妇血浆" && $config{$s}{"sample_info"}->[0] ne "地贫标准品");
	my $file = "$dOut/$s.genotype.txt";
	open(I, $file) or die $!;
	<I>;
	while(<I>){
		chomp;
		next if(/rs/);
		my @unit = split /\t/, $_;
		my ($cffdna, $genoM, $genoF, $p, $depa, $deph)=@unit[5,6,7,9,11,12];
		if(scalar @unit==17 || ($genoM ne "N/N" || $genoF ne "N/N")){
			## filter check
			my $is_filter = 0;
			my @finfo;
			if($cffdna < 0.05){
				$is_filter = 1;
				push @finfo, "LowCffdna";
			}
			if($depa < 1000){
				$is_filter = 1;
				push @finfo, "LowDep";
			}
			if($p < 0.75){
				$is_filter = 1;
				push @finfo, "Lowp";
			}
			my $check = scalar @finfo>0? join(",",@finfo): "Normal";
			if($#unit>13){
				print O join("\t",@unit[0..9], @unit[11,12], $check, @unit[13..$#unit]),"\n";
			}else{
				print O join("\t",@unit[0..9],@unit[11,12], $check),"\n";
			}
		}
	}
	close(I);
}
close(O);

# gather sample cffdna result
open (O,">$dOut/cffdna.txt") or die $!;
print O join("\t","Sample","SNPNum","AverDepth","AverDepthHotspot","%cffdna_EM",),"\n";
foreach my $s (sort {$a cmp $b} keys %config) {
	next if($config{$s}{"sample_info"}->[0] ne "地贫孕妇血浆" && $config{$s}{"sample_info"}->[0] ne "地贫标准品");
	my $file = "$dOut/$s.genotype.txt";
	my ($cffdna, $nsnp, $dep, $dephot)=("NA","NA","NA","NA");
	if (-f $file) {
		my $line = `cat $file|head -2|tail -1`;
		chomp $line;
		($cffdna, $nsnp, $dep, $dephot)=(split /\t/, $line)[5,10..12];
	}

	print O join("\t",$s,$nsnp,$dep,$dephot,$cffdna),"\n";
}
close (O) ;


#

######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ==============================================================
# sub function
# ==============================================================


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <chouxian.ma\@genetalks.com> 
Discription:
	
Usage:
  Options:
  -i		<file>	input samplesheet file, forced
  -id		<dir>	dir of input file(xxx.hotspots.result),xxx/06.hotspots_result/, forced
  -k		<str>	key of output file, forced
  -od		<dir>	output dir, optional, default [$dOut]
  -d        <int>   min depth, default [$min_depth]
  -e                infer genotype by cffdna
  -h			Help

USAGE
	print $usage;
	exit;
}

