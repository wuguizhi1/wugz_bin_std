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
my ($time, $sample, $primer, $samplename, $outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"t:s"=>\$time,
				"s:s"=>\$sample,
				"n:s"=>\$samplename,
				"od:s"=>\$outdir,
				"p:s"=>\$primer,
				) or &USAGE;
&USAGE unless ($time and $sample);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my$bamfile = "/data/bioit/biodata/duyp/Project/Thalassemia/PAGB_TMA_$time/05.variant_detect/$sample.GATK_realign.bam";
my$primerA = `samtools view $bamfile|less |grep 'SEA'|less |perl -ne '{chomp; (\$p)=\$_=~/XP:Z:(\\S+)/;print \$p,\"\\n\";}'|less|sort|uniq -c |awk '{print \$2}'`;
chomp($primerA);my @primerAll = split/\n/,$primerA;
foreach my $primer1(@primerAll){
	my $filename = "$time\_$sample.$primer1";
	if( defined $samplename ){  $filename = $samplename.".".$filename;  }
	my$cmd = "samtools view $bamfile|grep '$primerA'|less >$outdir/$filename.sam";
	print $cmd; system $cmd;

	$cmd = "perl $Bin/SEA_intersize.pl -i $outdir/$filename.sam -k $filename -od $outdir >$outdir/$filename.log";
	print $cmd,"\n"; system $cmd;
	`less $outdir/$filename.size |awk '{if(\$3>350){print \$1}}'|xargs |perl -ne '{
		chomp;
		my \@id = split;
		my %hash;
		foreach(\@id){\$hash{\$_}=1;}
		open(I, "samtools view -h $bamfile|") or die \$!; 
		open(OUT,"|samtools view -bS - >$outdir/$sample.bam");
		while(<I>){
			my \@line = split;
			if( not exists \$hash{\$line[0]} ){
				print OUT "\$_";
			}
		}
		close I;
		close OUT;
	}'`;
	`samtools sort $outdir/$sample.bam $outdir/$sample`;
	`samtools index $outdir/$sample.bam`;
	`perl /data/bioit/biodata/duyp/bin/ctDNA_dev/ctDNA/variant/fusion_detect.pl -i $outdir/$sample.bam -k $sample -od $outdir/`;
	#`rm $outdir/*bam*`;
	`rm $outdir/*log $outdir/*sam`;
	#$cmd = "Rscript /data/bioit/biodata/zenghp/bin/tools/drawing/density/multiDensity.r --infile $outdir/$filename.size --outfile $outdir/$filename.size.png  --value.col 3 --group.col 2 --group.lab \"group lab\" --x.lab \"Distance to BreakPoint\" --y.lab \"Frequency\" --title.lab \"$filename\" --skip 1 ";
	#print $cmd,"\n"; system $cmd;

	#$cmd = "Rscript hist.r --infile $outdir/$filename.size --outfile $outdir/$filename.size.hist.png --value.col 3 --group.col 2 --group.lab \"group lab\" --x.lab \"Distance to BreakPoint\" --y.lab \"Frequency\" --title.lab \"$filename\" --skip 1 ";
	#print $cmd,"\n"; system $cmd;
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
  -t	time, forced
  -s	sample name, forced
  -p	primer id, forced
  -h		 Help

USAGE
	print $usage;
	exit;
}

