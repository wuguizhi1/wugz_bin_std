#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

#use utf8;
#binmode(STDOUT, ":utf8");
#binmode(STDERR, ":utf8");
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($infile, $indir,$fkey,$outdir);
my $mode = 'e';
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$infile,
				"s:s"=>\$samplesheet,
				"id:s"=>\$indir,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($flist and $indir and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

open(IN, "$infile") or die "no such file: $infile\n";
while(<IN>){
	chomp;
	my @line = split/\t/, $_; 
	my( $SampleID, $SampleName, $Sample, $Geno, $Total, $Mutation, $Ratio, $cffDNA, $Genotype, $Pvalue ) = @line[0,1,7..14];
	if( $SampleID ne "SampleID" ){
		
	}
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath
{		#获取指定目录或文件的决定路径
		my ($type,$input) = @_;

		my $return;
		if ($type eq 'dir')
		{
				my $pwd = `pwd`;
				chomp $pwd;
				chdir($input);
				$return = `pwd`;
				chomp $return;
				chdir($pwd);
		}
		elsif($type eq 'file')
		{
				my $pwd = `pwd`;
				chomp $pwd;

				my $dir=dirname($input);
				my $file=basename($input);
				chdir($dir);
				$return = `pwd`;
				chomp $return;
				$return .="\/".$file;
				chdir($pwd);
		}
		return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i  <file>   	Input samplesheet file, forced
  -id  <dir>  	dir of xxx.variants.EM.txt, forced
  -k  <str>	Key of output file, forced
  -od <dir>	Dir of output file, default ./
  -m            mode, e or a, e: 提取指定位点的突变数据，a:输出所有样本的所有突变信息 
  -h		 Help

USAGE
	print $usage;
	exit;
}

