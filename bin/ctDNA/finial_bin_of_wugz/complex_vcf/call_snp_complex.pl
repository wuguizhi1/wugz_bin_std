#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($bamfile, $outdir, $outfile);
my ($print_type, $mindepth);
GetOptions(
"help|?" =>\&USAGE,
"i:s"=>\$bamfile,
"o:s"=>\$outfile,
"mindepth:s"=>\$mindepth,
"d:s"=>\$outdir,
"p:s"=>\$print_type,
) or &USAGE;
&USAGE unless defined($bamfile and $outfile);



my $RefGenome = "/data/bioit/biodata/duyp/bin/hg19/hg19.fasta";
if(not defined $print_type){    $print_type=0;  }
if(not defined $mindepth){      $mindepth=3;    }

$outdir||="./";
`mkdir $outdir` unless (-d $outdir);
$outdir=&AbsolutePath("dir",$outdir);

my(%all_bam_single, %all_bam_double);
open(I, "samtools view $bamfile|") or die $!;
while(<I>){
	chomp; my$lines=$_;
	my@line=split/\t/,$lines;

	my($id, $flag, $chr, $pos, $cigar, $seqSeq) = @line[0..3,5,9];
	next if( $cigar eq "*" );
	my $len = len_leave($cigar, "MD");
	my $end = $pos+$len-1;
	my $MD = "NA";
	if( $lines =~ m/MD:Z:(.+?)\t/){  
		$MD = $1; 
	}else{
		my $seqRef = get_samtools_faidx($chr,$pos,$end);
		$MD = trans_MD($seqSeq, $seqRef, $chr,$pos,$end, $cigar);
	}

	if( exists $all_bam_single{$id} ){
		my $each_flag = $all_bam_single{$id};
		foreach my$flag_(keys %$each_flag){
			my($chr_, $pos_, $end_, $cigar_, $seqSeq_, $MD_) = @{$all_bam_single{$id}{$flag_}};
			my($chrD, $posD, $endD, $Ins) = company_scale($chr,$pos,$end,$cigar, $chr_, $pos_, $end_, $cigar_);
			if( defined $endD ){
				push @{$all_bam_double{$id}}, ($flag,$flag_,$chrD,$posD,$endD,$Ins);
			}
		}
	}
	push @{$all_bam_single{$id}{$flag}}, ($chr,$pos,$end,$cigar,$seqSeq,$MD);
}
close I;

open(OUT,">$outfile");
foreach my$id(keys %all_bam_single){
	my$each_bam_id = $all_bam_single{$id};
	my($flag1,$flag2,$chrD,$posD,$endD,$Ins);
	if( exists $all_bam_double{$id} ){
		($flag1,$flag2,$chrD,$posD,$endD,$Ins) = @{$all_bam_double{$id}};
		print "$id\t$flag1,$flag2,$chrD,$posD,$endD,$Ins\n";
	}
=pob
	foreach my$flag(sort keys %$each_bam_id){
		my($chr,$pos,$end,$cigar,$seqSeq,$MD) = @{$all_bam_single{$id}{$flag}};
		if( $MD eq "NA" ){
			my $seqRef = get_samtools_faidx($chr,$pos,$end);
		}
	}
=cut
}
close OUT;
####################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub USAGE{
	my $usage=<<"USAGE";
	Program: $0
	Version: $version
	Contact:Wu guizhi<guizhi.wu\@genetalks.com>

	Usage:
		Options:
		-i      <file>  Input bam file, forced
		-o      <file>  Output file, forced
		-d      <str>   Outdir, [./]
		-h              Help

USAGE
print $usage;
exit;
}
sub trans_MD{  # 根据测序结果和参考序列得到MD
	my($seqSeq, $seqRef, $chr,$pos,$end, $cigar)=@_;
	my $re;
	if( $seqSeq =~ /$seqRef/ ){
		$re = $end-$pos+1;
	}else{
		my ($match_n_,$match_str_) = &sep_match_info($cigar);
		my @match_str = @$match_str_;
		foreach my$index(0..$#match_str){
			if( $match_str[$index] eq "I" ){ 
				$seqSeq = substr($seqSeq, $match_n_->[$index]);
			}elsif( $match_str[$index] eq "D" ){ 
				$re .= "^".substr($seqRef, 0, $match_n_->[$index]);
				$seqRef = substr($seqRef, $match_n_->[$index]);
			}elsif( $match_str[$index] eq "S" ){
				$seqSeq = substr($seqSeq, $match_n_->[$index]);
			}elsif( $match_str[$index] eq "M" ){
				my $refM = substr($seqRef, 0, $match_n_->[$index]);
				my $seqM = substr($seqSeq, 0, $match_n_->[$index]);
				if( $refM eq $seqM ){
					$re .= $match_n_->[$index];
				}else{
					my@refMa = split//,$refM;
					my@seqMa = split//,$seqM;
					my$len_each = 0;
					foreach my$m(0..$#refMa){
						if( $refMa[$m] eq $seqMa[$m] ){ $len_each++;
						}else{
							$re .= $len_each; $len_each = 0;
							$re .= $refMa[$m];
						}
					}
					if( $len_each > 0 ){ $re .= $len_each; }
				}
				$seqRef = substr($seqRef, $match_n_->[$index]);
				$seqSeq = substr($seqSeq, $match_n_->[$index]);
			}
		}
	}
	return($re);
}
sub company_scale{
	my($chr,$pos,$end,$cigar, $chr_, $pos_, $end_, $cigar_) = @_;
	if( $chr eq $chr_ ){
		if( $pos <= $pos_ and $end >=$pos_ ){
			if( $cigar =~ /I$/ ){
				return($chr, $pos_, $end, "endI");
			}elsif( $cigar_ =~ /^\d+I/ ){
				return($chr, $pos_, $end, "startI");
			}else{
				return($chr, $pos_, $end, "noI");
			}
		}elsif( $pos <= $end_ and $end >=$end_ ){
			if( $cigar_ =~ /I$/ ){
				return($chr, $pos, $end_, "endI");
			}elsif( $cigar =~ /^\d+I/ ){
				return($chr, $pos_, $end, "startI");
			}else{
				return($chr, $pos_, $end, "noI");
			}
		}
	}
	return(0);
}
sub get_samtools_faidx{  # 指定位置，截取基因组序列
	my($chr,$sta,$end)=@_;
	my $genome_ref = `samtools faidx $RefGenome $chr:$sta-$end`;
	$genome_ref =~ tr/ATGCatgc/ATGCATGC/;
	my ($id,@seq_ref) = split( /\n/,$genome_ref );
	if( scalar(@seq_ref) == 0 ){ print "what:$genome_ref"; }
	my $ref_seq12 = join "",@seq_ref;
	return($ref_seq12);
}
sub sep_match_info{   #分割cigar
	my($match_info)=@_;
	my @ucigar = split //, $match_info;
	my (@match_n,@match_str);
	my $I_before = 0;
	my $cigar="";
	foreach my $i(0..$#ucigar){
		if($ucigar[$i] eq "M" || $ucigar[$i] eq "I"|| $ucigar[$i] eq "D" || $ucigar[$i] eq "S"){
			push @match_str, $ucigar[$i];
			push @match_n, $cigar+$I_before;
			$cigar = "";
		}else{
			$cigar .= $ucigar[$i];
		}
	}
	return(\@match_n,\@match_str)
}
sub len_leave{  # 计数各个类型的长度
	my($cigar, $typ)=@_;
	my ($match_n_,$match_str_) = &sep_match_info($cigar);
	my @match_n = @$match_n_;
	my @match_str = @$match_str_;
	my ($reRef,$reSeq) = (0,0);
	foreach my$index(0..$#match_n){
		if( $match_str[$index] =~ /[MSI]/ ){
			$reSeq += $match_n[$index];
		}
		if( $match_str[$index] =~ /[MD]/ ){
			$reRef += $match_n[$index];
		}
	}
	if( $typ eq "MSI" ){
		return($reSeq);
	}elsif( $typ eq "MD" ){
		return($reRef);
	}
}
sub AbsolutePath{  #获取指定目录或文件的决定路径
	my ($type,$input) = @_;
	my $return;
	if ($type eq 'dir'){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	}elsif($type eq 'file'){
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
sub SHOW_TIME {  #显示当时时间函数，参数内容是时间前得提示信息，为字符串
	my ($str)=@_;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	my $temp=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
	print "$str:\t[".$temp."]\n";
}

