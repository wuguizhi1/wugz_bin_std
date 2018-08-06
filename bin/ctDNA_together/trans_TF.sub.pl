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

#######################################################################
##########  1：母亲 	2：胎儿		3：父亲		#######################
#######################################################################
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
#######################################################################

#######################################################################
#######################################################################
my(%Ratio_low, %Ratio_high, %Sample_Fetal_Geno);
open(IN, "$infile") or die "no such file: $infile\n";
open(OUT,">$outfile");
while(<IN>){
	chomp;
	my @line = split/\t/, $_; 
	my( $SampleID, $SampleName, $Sample, $Geno, $Total, $Mutation, $Ratio, $cffDNA ,$Genotype, $Pvalue ) = @line[0,1,6..13];
	$Geno = normal_Geno($Geno);
	my( $Mother, $Fetal ) = split/\|/, $Genotype;
	if( not defined $Fetal ){	next;	}

	my( $sample, $id ) = split/-/, $Sample;
	if( exists $people_geno{$sample}{"2"} ){
		my $Fetal_E = $people_geno{$sample}{"2"};
		my @Fetal_Geno = split/[\/,]/, $Fetal_E;
		foreach my$Fetal_Geno_each(@Fetal_Geno){
			$Sample_Fetal_Geno{$SampleName}{$Fetal_Geno_each} = 1;
		}
	}
	
	my @print_line = ($SampleName, $Sample, $Geno, $Total, $Mutation, $Ratio, $cffDNA , $Genotype, $Pvalue, $Mother, $Fetal);
	if( $SampleID ne "SampleID" and $Ratio > 0.01 ){
		$Sample_Fetal_Geno{$SampleName}{$Geno} = 1;
		if( not exists $Ratio_high{$SampleName}{$Geno} ){
			@{$Ratio_high{$SampleName}{$Geno}} = @print_line;
		}else{
			my @print_line_old = @{$Ratio_high{$SampleName}{$Geno}};
			if( $print_line[5] > $print_line_old[5] ){ 
				@{$Ratio_high{$SampleName}{$Geno}} = @print_line;
			}
		}
	}else{ 
		@{$Ratio_low{$SampleName}{$Geno}} = @print_line;
	}
}
close IN;
#######################################################################


foreach my$S(keys %Sample_Fetal_Geno){ 					#print "$S\n";
	my $Fetal_Geno = $Sample_Fetal_Geno{$S};
	foreach my$G( keys %$Fetal_Geno ){					#print "\t$G\n";
		if( exists $Ratio_high{$S}{$G} and $G ne "N" ){	
			my @print_line = @{$Ratio_high{$S}{$G}};	#print "\t$G\thigh\t@print_line\n";
			my( $re, $Mother_E, $Fetal_E, $Father_E ) = &Geno_JianYan(@print_line[0..2,9,10]);
			my $p = join "\t", (@print_line, $Mother_E, $Fetal_E, $Father_E);
			print OUT "$p\t$re\n";
		}elsif( exists $Ratio_low{$S}{$G} and $G ne "N" ){ 
			my @print_line = @{$Ratio_low{$S}{$G}};		#print "\t$G\tlow@print_line\n";
			my( $re, $Mother_E, $Fetal_E, $Father_E ) = &Geno_JianYan(@print_line[0..2,9,10]);
			my $p = join "\t", (@print_line, $Mother_E, $Fetal_E, $Father_E);
			print OUT "$p\tFN\n";
		}
	}
}
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
sub Geno_JianYan{
	my( $SampleName, $Sample, $Geno, $Mother, $Fetal ) = @_;
	my( $sample, $id ) = split/-/, $Sample;
	#print "wugz:$SampleName, $sample, $Geno, $Mother, $Fetal\n";
	if( exists $people_geno{$sample} ){
		my( $re, $Mother_E, $Fetal_E, $Father_E ) = ("NA", "N/N", "N/N", "N/N");
		if( exists $people_geno{$sample}{"1"} ){	
			$Mother_E = $people_geno{$sample}{"1"};
		}
		if( exists $people_geno{$sample}{"2"} ){	
			$Fetal_E = $people_geno{$sample}{"2"};
		}
		if( exists $people_geno{$sample}{"3"} ){	
			$Father_E = $people_geno{$sample}{"3"};
		}
		my $Geno_E = $Mother_E.$Fetal_E;
		#print "wugz:$re, $Mother_E, $Fetal_E, $Father_E\n";

		if( $Geno_E =~ /\Q$Geno\E/ ){	
			if( $Fetal =~ /\Q$Geno\E/ ){ # 胎儿阳性
				if( $Fetal_E =~ /$Fetal/ ){
					$re = "TP";
				}else{
					$re = "FP";
				}
			}else{	# 胎儿阴性
				if( $Fetal_E =~ /\Q$Geno\E/ ){
					$re = "FN";
				}else{
					$re = "TN";
				}
			}
		}else{
			$re = "NA";
		}
		#print "$Mother_E.$Fetal_E\tP:(Fetal)$Fetal###(Geno)$Geno\tT:$Fetal##(Fetal_E)$Fetal_E\t$re\n";
		return($re, $Mother_E, $Fetal_E, $Father_E);
	}else{
		print "no sample named $sample in $samplelist!\n";
	}
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

