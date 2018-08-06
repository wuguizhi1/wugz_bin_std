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
my ($infile, $samplelist, $outfile);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$infile,
				"s:s"=>\$samplelist,
				"od:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($infile and $samplelist and $outfile);

my %people_geno;
open(PEOPLE, $samplelist) or die "no such file: $samplelist!\n";
while(<PEOPLE>){
	chomp;
	my@line = split/\t/, $_;
	my( $sample, $id ) = split/-/, $line[1];
	if( not defined $id ){	next;	}
	if( not defined $id ){	print "@line\n";next;	}
	if( $line[8] ne "" and scalar(@line) == 9 ){
		 $people_geno{$sample}{$id} = $line[8];
	}elsif( $line[9] ne "" and scalar(@line) == 10 ){
		$people_geno{$sample}{$id} = $line[9];
	}else{
		#print "#$line[0]\t#$line[1]\t#$line[2]\t#$line[3]\n";
	}
}
close PEOPLE;

my(%Ratio_low, %Ratio_high, %Sample_Fetal_Geno, $print_pre);
open(IN, "$infile") or die "no such file: $infile\n";
while(<IN>){
	chomp; my$lines = $_;
	my @line = split/\t/, $lines; 
	@line = @line[0..17]; $line[8] =~ s/^'-/-/g; $line[14] =~ s/^'-/-/g; $line[15] =~ s/^'-/-/g;
	my($SampleName, $Sample, $MutID, $Ratio, $MaternalGenotype, $FeatlGenotype ) = @line[0,1,8,11,14,15];
	if( $MutID =~ /rs/ ){
		$print_pre .= "$lines\n";
	}elsif( $MutID =~ /MutID/ ){
		$print_pre .= "$lines\tFetalGiving\tMaternalGiving\tPaternalGiving\tTrueOrFalse\n";
	}else{
		my( $sample, $id ) = split/-/, $Sample;
		my($Fetal, $Mother, $Father) = Geno_E($SampleName, $sample);
		@line = (@line, $Fetal, $Mother, $Father);
		my $Geno = normal_Geno($MutID); 
		$line[8] = $Geno;
		$line[14] = normal_Geno($MaternalGenotype);
		$line[15] = normal_Geno($FeatlGenotype);
		if( $Ratio > 0.01 ){
			$Sample_Fetal_Geno{$SampleName}{$Geno} = 1;
			if( not exists $Ratio_high{$SampleName}{$Geno} ){
				@{$Ratio_high{$SampleName}{$Geno}} = @line;
			}else{
				my @line_old = @{$Ratio_high{$SampleName}{$Geno}};
				if( $line[11] > $line_old[11] ){
					@{$Ratio_high{$SampleName}{$Geno}} = @line;
				}
			}
		}else{
			@{$Ratio_low{$SampleName}{$Geno}} = @line;
		}
	}
}
close IN;

open(OUT,">$outfile");
print OUT $print_pre;
foreach my$S(keys %Sample_Fetal_Geno){
	my %print_G;
	my $Fetal_Geno = $Sample_Fetal_Geno{$S};
	foreach my$G( keys %$Fetal_Geno ){
		my @print_line;
		if( exists $Ratio_high{$S}{$G} and $G ne "N" ){
			@print_line = @{$Ratio_high{$S}{$G}};
			$print_G{$S}{$G} = 1;
		}elsif( exists $Ratio_low{$S}{$G} and $G ne "N" ){ 
			@print_line = @{$Ratio_low{$S}{$G}};
			$print_G{$S}{$G} = 1;
		}
		
		if( @print_line ){
			my $re = &Geno_JianYan(@print_line[8,14,15,18..20]);
			foreach my$i(0..$#print_line){$print_line[$i] =~ s/^-/'-/g;}
			my $p = join "\t", @print_line;
			print OUT "$p\t$re\n";
		}else{}
	}
	
	my $low_Geno = $Ratio_low{$S};
	foreach my$G(keys %$low_Geno){
		if( not exists $print_G{$G} ){
			my@print_line = @{$Ratio_low{$S}{$G}};
			foreach my$i(0..$#print_line){$print_line[$i] =~ s/^-/'-/g;}
			my $p = join "\t", @print_line;
			print OUT "$p\n";
		}
	}
}
close OUT;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub normal_Geno{
	my( $Geno ) = @_;
	$Geno =~ s/-28_C/-28/g;
	$Geno =~ s/-28_G/-28/g;
	$Geno =~ s/CD71-72_A/CD71-72/g;
	$Geno =~ s/CD71-72_T/CD71-72/g;
	return($Geno);
}
sub Geno_E{
	my( $SampleName, $sample ) = @_;
	my( $Mother_E, $Fetal_E, $Father_E ) = ("N/N", "N/N", "N/N");
	if( exists $people_geno{$sample} ){
		my @Fetal_Geno ;
		if( exists $people_geno{$sample}{"1"} ){
			$Mother_E = $people_geno{$sample}{"1"};
		}
		if( exists $people_geno{$sample}{"2"} ){
			$Fetal_E = $people_geno{$sample}{"2"};
			@Fetal_Geno = split/[\/,]/, $Fetal_E;
			foreach my$Fetal_Geno_each(@Fetal_Geno){
				$Sample_Fetal_Geno{$SampleName}{$Fetal_Geno_each} = 1;
			}
		}
		if( exists $people_geno{$sample}{"3"} ){
			$Father_E = $people_geno{$sample}{"3"};
		}
	}else{
		print "no sample named $sample in $samplelist!\n";
	}
	return( $Fetal_E, $Mother_E, $Father_E);
}
sub Geno_JianYan{
	my($Geno, $Mother, $Fetal, $Fetal_E, $Mother_E, $Father_E, $re) = @_; 
	my $Geno_E = $Mother_E.$Fetal_E;

	if( $Geno_E =~ /\Q$Geno\E/ ){    # 测得基因型来源（母亲与胎儿）	
		if( $Fetal =~ /\Q$Geno\E/ ){ # 胎儿阳性
			if( $Fetal_E =~ /\Q$Fetal\E/ and $Mother_E =~ /\Q$Mother\E/ ){
				$re = "TP";
			}elsif( $Fetal_E =~ /\Q$Fetal\E/ and $Mother eq "N/N" and $Mother_E !~ /\Q$Geno\E/){
				$re = "TP";
			}else{
				$re = "FP";
			}
		}else{	# 胎儿阴性
			if( $Fetal_E =~ /\Q$Geno\E/  ){
				$re = "FN";
			}elsif( $Mother_E =~ /\Q$Geno\E/ and $Mother_E =~ /\Q$Mother\E/ ){
				$re = "TN";
			}elsif( $Mother_E !~ /\Q$Geno\E/ and $Mother !~ /\Q$Geno\E/ ){
				$re = "TN";
			}else{
				$re = "FN";
			}
		}
	}else{	#得基因型来源（父亲），或未给定
		$re = "FP(?)";
	}
	print "$re:\tM:$Mother\tf:$Fetal\t\t\tM:$Mother_E\tf:$Fetal_E\tF:$Father_E\t\t$Geno\n";
	return($re);
}

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

