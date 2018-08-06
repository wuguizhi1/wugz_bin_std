#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($samplelist, $inputdir, $cffDNA, $fkey, $outdir);
my $max_x = 400;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$samplelist,
				"id:s"=>\$inputdir,
				"c:s"=>\$cffDNA,
				"m:s"=>\$max_x,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($samplelist and $fkey and $cffDNA and $inputdir);
$outdir||="./";
mkdir $outdir if (! -d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my @bam_info;
open(IN, $samplelist) or die "no such file: $samplelist!\n";
while(<IN>){
	chomp;
	push @bam_info,$_;
}
close IN;
my $fsize="$outdir/$fkey.size.txt";
`rm $fsize` if(-f $fsize);
foreach my $s(@bam_info){
	my $bam = "$inputdir/$s.bam";  
	my $cmd ="samtools view -f 0x40 -F 0x800 $bam | perl -ane '{if(\$F[8]!=0 && abs(\$F[8])<$max_x){print \"$s\",\"\\t\",abs(\$F[8]),\"\\n\"}}' >>$fsize";
	print "$cmd\n";
	`$cmd`;
}
my $cmdR = "Rscript $Bin/InsertSize.r 200 400 $cffDNA $inputdir $fsize $fsize.adjust.txt $outdir";
print $cmdR,"\n";
`$cmdR`;

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
	$/="\n";

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
Contact: zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i     <str,>   bam infos, "sampleName", forced
  -id	 <str>    dir of this samplesheet, forced
  -c     <str>    file of cffDNA, forced
  -m     <int>    max x, [400]
  -k     <str>    key of output file, forced
  -od    <dir>    output dir, optional

  -h                     Help

USAGE
	print $usage;
	exit;
}

