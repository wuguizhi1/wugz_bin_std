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
my ($flist, $indir,$fkey,$outdir);
my $mode = 'e';
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$flist,
				"id:s"=>\$indir,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				"m:s"=>\$mode,
				) or &USAGE;
&USAGE unless ($flist and $indir and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my %info;
my @f = glob("$indir/snp_indel/*.variants.EM.txt");
foreach my $f(@f){
	my ($s)=$f=~/$indir\/snp_indel\/(\S+).variants.EM.txt/;
	open(F,$f) or die $!;
	<F>;
	while(<F>){
		chomp;
		my @unit = split;
		my ($geno, $t, $m, $r, $f, $g, $p)=@unit[5..11];
		next if ($geno=~/^rs/);
		$geno=~s/\w+://;
		@{$info{$s}{$geno}}=($t,$m,$r,$f,$g, $p);
		# Total Mutation Ratio fetal.EM Refined.genotype.EM p
	}
	close(F);
}
#SampleName      ID      N       b       f
#L1803300331     SEA     1554    0       0.000000
#L1803300332     SEA     1749    0       0.000000
open(S, "$indir/alpha.SEA_fusion.result.xls") or die $!;
<S>;
while(<S>){
	chomp;
	my ($s, $g, $t, $m, $r, $f, $gt, undef, $p) = split /\t+/, $_;
	@{$info{$s}{$g}}=($t,$m,$r, $f, $gt, $p); 
	# TotalDepth  FusionDepth  FusionFraction  cffdna_SNP  MaxProbGenotype MaxProb
}
close(S);

my %geno_name=("--SEA"=>"SEA",
		"αα"=>"SEA",
		"IVS-II-5"=>"+5",
		"CD27-28"=>"CD27/28",
		"CD71-72"=>"CD71-72_A",
		"-28"=>"-28_G",
		"βE" => 'CD26',
		"αWS" => "WS",		
);


open(L,$flist) or die $!;
#binmode(L, ':encoding(utf8)');
open(O,">$outdir/$fkey.result") or die $!;
my $title = <L>;
chomp $title;
print O join("\t",$title,"Geno","Total","Mutation","Ratio", "cffDNA", "Genotype", "Pvalue"),"\n";
if($mode eq 'e'){
while(<L>){
	chomp;
	next if(/^$/);
	my @unit = split /\t/, $_;
	my ($s, $geno)=@unit[1,7];
	if(!defined $geno){
		print $_,"\n";
		print "no geno\n";
		die;
	}
	$geno=~s/^,//;
	$geno=~s/--//;
	my @geno =split /,/, $geno;
	foreach my $g (@geno){
		next if ($g =~/αα\/αα/);
		next if ($g =~/β.?N/ || $g =~/α.?N/);
		($g)=split /\//, $g if (!exists $geno_name{$g});
		#$g=~s/\/N//;
		my $g_new = exists $geno_name{$g}?  $geno_name{$g}: $g;
		if(!exists $info{$s}{$g_new}){
			#print Dumper %{$info{$s}};
			print $_,"\n";
			print $s,"\t", $g, "\t", $g_new, "\n";
			print O join("\t", $_, $g, "NA", "NA", "NA", "NA", "NA", "NA"),"\n";
			next;
		}
		print O join("\t", $_, $g, @{$info{$s}{$g_new}}),"\n";
	}
}
}elsif($mode eq 'a'){
while (<L>){
	chomp;
	next if (/^$/);
	my @unit = split /\t/, $_;
	my ($s, $geno)=@unit[1,7];
	
	foreach my $g (keys %{$info{$s}}){
		print O join("\t", $_, $g, @{$info{$s}{$g}}),"\n";
	}
}
}else{
	die "mode error!!";
}
close(L);
close(O);


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

