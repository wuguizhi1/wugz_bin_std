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
my ($infile, $sampleGeno, $outfile);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$infile,
				"g:s"=>\$sampleGeno,
				"o:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($infile and $sampleGeno and $outfile);

#######################################################################
my %Sample_Fetal_Geno;
my($Fetal, $Mother, $Father) = split":", $sampleGeno;   print "f:$Fetal M:$Mother F:$Father\n"; 
my @Fetal_Geno = split/[\/,]/, $Fetal;
my @Mother_Geno = split/[\/,]/, $Mother;
foreach my$Fetal_Geno_each(@Fetal_Geno,@Mother_Geno){
	$Sample_Fetal_Geno{$Fetal_Geno_each} = 1;
}
my(%Ratio_low, %Ratio_high);
open(IN, "$infile") or die "no such file: $infile\n";
my $print_pre;
while(<IN>){
	chomp;
	my $lines = $_;
	my @line = split/\t/, $lines; 
	@line = @line[0..12];
	$line[1] =~ s/^'-/-/g;
	$line[6] =~ s/^'-/-/g;
	$line[7] =~ s/^'-/-/g;
	my( $Sample, $MutID, $Ratio, $MaternalGenotype, $FeatlGenotype ) = @line[0..1,4,6..7];
	if( $MutID =~ /rs/ ){
		$print_pre .= "$lines\n";
	}elsif( $MutID =~ /MutID/ ){
		$print_pre .= "$lines\tFetalGiving\tMaternalGiving\tPaternalGiving\tTrueOrFalse\n"; 
	}else{
		@line = (@line, $Fetal, $Mother, $Father);
		my $Geno = normal_Geno($MutID);
		$line[1] = $Geno;
		$line[6] = normal_Geno($MaternalGenotype);
		$line[7] = normal_Geno($FeatlGenotype);
		if( $Ratio > 0.01 ){
			$Sample_Fetal_Geno{$Geno} = 1;
			if( not exists $Ratio_high{$Geno} ){
				@{$Ratio_high{$Geno}} = @line;
			}else{
				my @line_old = @{$Ratio_high{$Geno}};
				if( $line[4] > $line_old[4] ){
					@{$Ratio_high{$Geno}} = @line;
				}
			}
		}else{
			@{$Ratio_low{$Geno}} = @line;
		}
	}
}
close IN;
#######################################################################
open(OUT,">$outfile");
print OUT $print_pre;
my %print_G;
foreach my$G(keys %Sample_Fetal_Geno){ 					
	my @print_line;
	if( exists $Ratio_high{$G} and $G ne "N" ){	
		@print_line = @{$Ratio_high{$G}};	
		$print_G{$G} = 1;
		#print "\t$G\thigh\t@print_line\n";
	}elsif( exists $Ratio_low{$G} and $G ne "N" ){ 
		@print_line = @{$Ratio_low{$G}};
		$print_G{$G} = 1;	
		#print "\t$G\tlow@print_line\n";
	}

	if( @print_line ){
		my $re = &Geno_JianYan(@print_line[1,6..7,13..15]);
		foreach my$i(0..$#print_line){		$print_line[$i] =~ s/^-/'-/g;		}
		my $p = join "\t", @print_line;
		print OUT "$p\t$re\n";
		#print "$p\n";
	}else{
	}
}
foreach my$G(keys %Ratio_low){ 
	if( not exists $print_G{$G} ){
		my@print_line = @{$Ratio_low{$G}};
		foreach my$i(0..$#print_line){		$print_line[$i] =~ s/^-/'-/g;		}
		my $p = join "\t", @print_line;
		print OUT "$p\n";
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
  -i  <file>	Input sample file, forced
  -g  <str>	Iutput sampleGenotype, forced
  -o  <file>	Output file , forced
  -h		Help

USAGE
	print $usage;
	exit;
}

