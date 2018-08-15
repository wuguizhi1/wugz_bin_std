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


#################### read bamfile
my(%all_bam_single, %all_bam_double, %genome_ref, %genome_alt, %reDP, %reAD, %reADR);
open(I, "samtools view $bamfile|") or die $!;
while(<I>){
	chomp; my$lines=$_;
	my@line=split/\t/,$lines;

	my($id, $flag, $chr, $pos, $cigar, $seqSeq) = @line[0..3,5,9];
	next if( $cigar eq "*" );
	my $len = len_leave($cigar, "MD");
	my $end = $pos+$len-1;

	my $MD;
	if( $lines =~ m/MD:Z:(.+?)\t/){  
		$MD = $1; 
	}else{
		my $seqRef = get_samtools_faidx($chr,$pos,$end);
		$MD = trans_MD($seqSeq, $seqRef, $chr,$pos,$end, $cigar);
	}
	
	my $primer;
	if( $lines =~ m/XP:Z:(.+?)\t/){
		$primer = $1;
	}else{
		print "Error: no primer info! $lines\n"; next; 
	}

	next if( not defined $MD or not defined $primer );
	my($rnum, $is_reverse ,$is_unmap)=&explain_bam_flag($flag);
	&rember_genome($chr,$pos,$MD,$id);
	&counts_genome($id,$chr,$pos,$cigar,$seqSeq,$MD,$primer,$is_reverse);

	if( $cigar =~ /^(\d+)S/ or $cigar =~ /(\d+)S$/ ){
		my ($match_n_,$match_str_) = &sep_match_info($cigar);
		my @match_str = @$match_str_;
		my( $startN, $startL , $endN, $endL ) = (0, 0, $#match_str, 0);
		if( $match_str[0] eq "S" ){   $startN = 1; $startL=$match_n_->[0];  }
		if( $match_str[$#match_str] eq "S" ){ $endN = $#match_str-1; $endL=$match_n_->[$#match_str]; }
		$cigar = "";
		foreach my$i( $startN..$endN ){ $cigar .= $match_n_->[$i].$match_str[$i]; }
		$seqSeq = substr($seqSeq, $startL, length($seqSeq)-$startL-$endL);
	}

	if( exists $all_bam_single{$id} ){
		my $each_flag = $all_bam_single{$id};
		foreach my$flag_(keys %$each_flag){
			my($chr_, $pos_, $end_, $cigar_, $seqSeq_, $MD_, $primer_, $is_reverse_) = @{$all_bam_single{$id}{$flag_}};
			my($chrD,$posD,$endD,$cigar1,$cigar2,$seqSeqD1,$seqSeqD2) = 
            company_scale($chr,$pos,$end,$cigar,$seqSeq,$chr_,$pos_,$end_,$cigar_,$seqSeq_);
			if( defined $endD and $primer_ eq $primer ){
				my $typ;
				if( $cigar1 eq $cigar2 and $seqSeqD1 eq $seqSeqD2 ){ $typ = "doublesame";
				}else{  $typ = "doublecigarNe";
					my($cigarL, $cigarR, $seqSeqDL, $seqSeqDR);
					my %seq_double_hash;
					my $seq_double = find_diff_cigar($chrD,$posD,$endD,$cigar1,$cigar2,$seqSeqD1,$seqSeqD2,\%seq_double_hash);
					%seq_double_hash = %$seq_double;
					foreach my$key(keys %seq_double_hash){
						print "\tnot\t$key\t$seq_double_hash{$key}\n";
					}
				
				print $typ,join"\t",($chrD,$posD,$endD,$cigar1,$cigar2,$seqSeqD1,$seqSeqD2), "\n";
				print "\t", join"\t",($chr,$pos,$end,$cigar,$seqSeq), "\n";
				print "\t", join"\t",($chr_,$pos_,$end_,$cigar_,$seqSeq_), "\n";
				}
				push @{$all_bam_double{$id}}, ($typ, $flag, $flag_, $chrD, $posD, $endD);
			}
		}
	}
	push @{$all_bam_single{$id}{$flag}}, ($chr,$pos,$end,$cigar,$seqSeq,$MD,$primer,$is_reverse);
}
close I;


#################  count vcf
open(OUT,">$outfile");
foreach my$id(keys %all_bam_double){
	my($typ,$flag1,$flag2,$chrD,$posD,$endD,$Ins);
	if( exists $all_bam_double{$id} ){
		($typ,$flag1,$flag2,$chrD,$posD,$endD,$Ins) = @{$all_bam_double{$id}};
		my($chr1, $pos1, $end1, $cigar1, $seqSeq1, $MD1, $primer1, $is_reverse1) = @{$all_bam_single{$id}{$flag1}};
		my($chr2, $pos2, $end2, $cigar2, $seqSeq2, $MD2, $primer2, $is_reverse2) = @{$all_bam_single{$id}{$flag2}};
		print "$id\t$flag1,$flag2,$chrD,$posD,$endD\n";
	}
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################


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

sub rember_genome{
	my($chr,$pos,$MD,$id)=@_;
	my@md = split//, $MD;
	my($md_n, $md_str, $md_d, $md_D);
	foreach my$index(0..$#md){
		if( $md[$index] =~ /\d/ ){ 
			if( $md_D == 1 ){ $genome_alt{$chr}{$pos-1}{$id."D"} = $md_d; }
			$md_n .= $md[$index]; $md_D = 0;
		}elsif( $md[$index] =~ /\^/ ){
			$pos += $md_n; $md_n = 0; $md_D = 1
		}else{
			$pos += $md_n+1; $md_n = 0;
			$genome_ref{$chr}{$pos-1} = $md[$index];
			if( $md_D == 1 ){ $md_d .= $md[$index]; } 
			if( $md_D == 0 ){ $genome_alt{$chr}{$pos-1}{$id."M"} = $md[$index]; }
			#print "genome_ref: $chr\t$pos-1\t$md[$index]\n";
		}
	}
}
sub counts_genome{
	my($id,$chr,$pos,$cigar,$seqSeq,$MD,$primer,$is_reverse)=@_;
	my($pos_r, $pos_l) = (0, 0);
	my ($match_n_,$match_str_) = &sep_match_info($cigar);
	my @match_str = @$match_str_;
	foreach my$index(0..$#match_str){
		if( $match_str[$index] eq "M" ){
			my @seqM = split//, substr($seqSeq, $pos_l, $match_n_->[$index]);
			foreach my$i(0..($match_n_->[$index]-1)){
				$reDP{$chr}{$pos+$pos_r+$i}{$primer}{$seqM[$i]}++;
				$reAD{$chr}{$pos+$pos_r+$i}{$primer}{$seqM[$i]}++;
				$reADR{$chr}{$pos+$pos_r+$i}{$primer}{$seqM[$i]} += $is_reverse;
				if( exists $genome_alt{$chr}{$pos+$pos_r+$i}{$id."M"} ){
					$genome_alt{$chr}{$pos+$pos_r+$i}{$id."M"} = $seqM[$i];
				}
			}
			$pos_r += $match_n_->[$index]; $pos_l += $match_n_->[$index];
		}elsif( $match_str[$index] eq "D" ){
			my $Del = "";
			foreach my$i(0..($match_n_->[$index]-1)){ 
				if( not exists $genome_ref{$chr}{$pos+$pos_r+$i} ){
					$genome_ref{$chr}{$pos+$pos_r+$i} = get_samtools_faidx($chr,$pos+$pos_r+$i,$pos+$pos_r+$i);
				}
				$Del .= $genome_ref{$chr}{$pos+$pos_r+$i}; 
			}
			$reDP{$chr}{$pos+$pos_r-1}{$primer}{"D$Del"}++;
			$reAD{$chr}{$pos+$pos_r-1}{$primer}{"D$Del"}++;
			$reADR{$chr}{$pos+$pos_r-1}{$primer}{"D$Del"} += $is_reverse;
			$pos_r += $match_n_->[$index];
		}elsif( $match_str[$index] eq "I" ){
			my $Ins = substr($seqSeq, $pos_l, $match_n_->[$index]);
			$reDP{$chr}{$pos+$pos_r}{$primer}{"I$Ins"}++;
			$reAD{$chr}{$pos+$pos_r}{$primer}{"I$Ins"}++;
			$reADR{$chr}{$pos+$pos_r}{$primer}{"I$Ins"} +=$is_reverse;
			$genome_alt{$chr}{$pos+$pos_r}{$id."I"} = $Ins;
			$pos_l += $match_n_->[$index];
		}elsif( $match_str[$index] eq "S" ){
			$pos_l += $match_n_->[$index];
		}
	}
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
	my($chr,$pos,$end,$cigar,$seqSeq, $chr_, $pos_, $end_, $cigar_, $seqSeq_) = @_;
	my($posD,$endD,$cigar1,$cigar2,$seqSeqD1,$seqSeqD2);
	if( $chr eq $chr_ ){
		if( $pos < $pos_ and $end <=$end_ and $pos_ <= $end){
			# |----------------------|
			#        |----------------------|
			my ($cigar_1,$cigar_2) = &sep_cigar($cigar, $pos_-1, $pos); 
			my ($cigar_3,$cigar_4) = &sep_cigar($cigar_, $end, $pos_);
			my $len = len_leave($cigar_1, "MSI"); 
			my $len_ = len_leave($cigar_3, "MSI"); 
			$seqSeqD1 = substr($seqSeq, $len);
			$seqSeqD2 = substr($seqSeq_, 0, $len_);
			$posD = $pos_; $endD = $end; $cigar1 = $cigar_2; $cigar2 = $cigar_3; 
		}elsif( $pos > $pos_ and $end >=$end_ and $pos <= $end_ ){
			#        |----------------------|
			# |----------------------|
			my ($cigar_1,$cigar_2) = &sep_cigar($cigar, $end_, $pos);
			my ($cigar_3,$cigar_4) = &sep_cigar($cigar_, $pos-1, $pos_);
			my $len = len_leave($cigar_1, "MSI");
			my $len_ = len_leave($cigar_3, "MSI");
			$seqSeqD1 = substr($seqSeq, 0, $len);
			$seqSeqD2 = substr($seqSeq_, $len_);
			$posD = $pos; $endD = $end_; $cigar1 = $cigar_1; $cigar2 = $cigar_4;
		}elsif( $pos <= $pos_ and $end >=$end_  ){
			# |------------------------------------------|
			#  |----------------------|
			my ($cigar_1,$cigar_2) = &sep_cigar($cigar, $pos_-1, $pos);
			if( $pos == $pos_ ){ $cigar_2 = $cigar_1; }
			my ($cigar_3,$cigar_4) = &sep_cigar($cigar_2, $end_, $pos_);
			my $len = len_leave($cigar_1, "MSI");
			my $len_ = len_leave($cigar_3, "MSI");
			$seqSeqD1 = substr($seqSeq, $len, $len_);
			$seqSeqD2 = $seqSeq_; 
			if( $pos == $pos_ ){ $seqSeqD1 = substr($seqSeq, 0, $len_); }
			$posD = $pos_; $endD = $end_; $cigar1 = $cigar_3; $cigar2 = $cigar_;
		}elsif( $pos >= $pos_ and $end <$end_ ){
			#  |----------------------|
			# |------------------------------------------|
			my ($cigar_1,$cigar_2) = &sep_cigar($cigar_, $pos-1, $pos_);
			if( $pos == $pos_ ){ $cigar_2 = $cigar_1; }
			my ($cigar_3,$cigar_4) = &sep_cigar($cigar_2, $end, $pos);
			my $len = len_leave($cigar_1, "MSI");
			my $len_ = len_leave($cigar_3, "MSI");
			$seqSeqD1 = $seqSeq;
			$seqSeqD2 = substr($seqSeq_, $len, $len_);
			if( $pos == $pos_ ){ $seqSeqD2 = substr($seqSeq_, 0, $len_); }
			$posD = $pos; $endD = $end; $cigar1 = $cigar; $cigar2 = $cigar_3;
		}else{
			return(0);
		}
		return($chr,$posD,$endD,$cigar1,$cigar2,$seqSeqD1,$seqSeqD2);
	}
	return(0);
}
sub find_diff_cigar{
	my($chrD,$posD,$endD,$cigar_a,$cigar_b, $seqSeqD1,$seqSeqD2,$seq_double_hash)=@_;
	if( $cigar_a ne $cigar_b ){
		my($match_n_a,$match_str_a) = &sep_match_info($cigar_a);
		my($match_n_b,$match_str_b) = &sep_match_info($cigar_b);
		my $scalar_a = scalar(@$match_n_a);my $scalar_b = scalar(@$match_n_b);
		my( $start, $end ) = (1, scalar(@$match_n_a));
		if( scalar(@$match_n_b) < scalar(@$match_n_a) ){ $end = scalar(@$match_n_b); }
		my($cigar_l, $cigar_r) = ("", "");
		for( my$index=0; $index++; $index<$end ){
			if( $match_n_a->[$index] == $match_n_b->[$index] and $match_str_a->[$index] eq $match_str_b->[$index] ){
				$cigar_l .= $match_n_a->[$index].$match_str_a->[$index];
			}elsif( $match_str_a->[$index] eq $match_str_b->[$index] and $match_n_a->[$index] < $match_n_b->[$index] ){
				$cigar_l .= $match_n_a->[$index].$match_str_a->[$index]; last;
			}elsif( $match_str_a->[$index] eq $match_str_b->[$index] and $match_n_a->[$index] > $match_n_b->[$index] ){
				$cigar_l .= $match_n_b->[$index].$match_str_a->[$index]; last;
			}else{}
		}
		for( my$index=0; $index++; $index<$end ){
			my$indexa = $scalar_a-$index-1;my$indexb = $scalar_b-$index-1;
			if( $match_n_a->[$indexa] == $match_n_b->[$indexb] and $match_str_a->[$indexa] eq $match_str_b->[$indexb] ){
				$cigar_r .= $match_n_a->[$indexa].$match_str_a->[$indexa];
			}elsif( $match_str_a->[$indexa] eq $match_str_b->[$indexb] and $match_n_a->[$indexa] < $match_n_b->[$indexb] ){
				$cigar_r .= $match_n_a->[$indexa].$match_str_a->[$indexa]; last;
			}elsif( $match_str_a->[$indexa] eq $match_str_b->[$indexb] and $match_n_a->[$indexa] > $match_n_b->[$indexb] ){
				$cigar_r .= $match_n_b->[$indexb].$match_str_a->[$indexa]; last;
			}else{}
		}
		my$seqSeqDL1 = substr($seqSeqD1, 0, len_leave($cigar_l, "MSI") );
		my$seqSeqDL2 = substr($seqSeqD2, 0, len_leave($cigar_l, "MSI") );
		if( $seqSeqDL1 ne $seqSeqDL2 ){ 
			$seq_double_hash = find_diff_seqSeq($chrD,$posD,$endD,$cigar_l, $seqSeqDL1, $seqSeqDL1, $seq_double_hash);
		}else{
			my $endNew = $posD + len_leave($cigar_l, "MD") - 1;
			$seq_double_hash->{$chrD.":".$posD.":".$endNew.":".$cigar_l.":".$seqSeqDL1} = 1; 
			$posD = $endD - len_leave($cigar_r, "MD") + 1;
		}
		my$seqSeqDR1 = substr($seqSeqD1, length($seqSeqD1) - len_leave($cigar_r, "MSI") );
		my$seqSeqDR2 = substr($seqSeqD2, length($seqSeqD2) - len_leave($cigar_r, "MSI") );
		if( $seqSeqDR1 eq $seqSeqDR2 and $seqSeqDR1 !~ /N/ ){
			$seq_double_hash->{$chrD.":".$posD.":".$endD.":".$cigar_r.":".$seqSeqDR1} = "seqSame0";
		}elsif( length($seqSeqD1) == 1 ){
		}else{
			$seq_double_hash = find_diff_seqSeq($chrD,$posD,$endD,$cigar_r, $seqSeqDR1, $seqSeqDR1, $seq_double_hash);
		}
	}elsif( $seqSeqD1 eq $seqSeqD2 and $seqSeqD1 !~ /N/ ){
		$seq_double_hash->{$chrD.":".$posD.":".$endD.":".$cigar_a.":".$seqSeqD1} = "seqSame0";
	}elsif( length($seqSeqD1) == 1 ){
	}else{
		$seq_double_hash = find_diff_seqSeq($chrD,$posD,$endD,$cigar_a, $seqSeqD1, $seqSeqD2, $seq_double_hash);
	}
	return($seq_double_hash);
}
sub find_diff_seqSeq{
	my($chrD,$posD,$endD, $cigar,$seqSeqD1,$seqSeqD2,$seq_double_hash ) = @_;
	my @seq1 = split//,$seqSeqD1;
	my @seq2 = split//,$seqSeqD2;
	my($seqL, $n) = ("", 1);
	foreach my$index( 0..$#seq1 ){
		if( $seq1[$index] ne $seq2[$index] ){ $n = $index; last;
		}elsif( $seq1[$index] eq "N" or $seq2[$index] eq "N" ){ $n = $index; last;
		}elsif( $seq1[$index] eq $seq2[$index] ){ $seqL .= $seq1[$index];
		}else{ }
	}
	my ($cigar_1,$cigar_2) = &sep_cigar($cigar, $n+$#seq1-1, $#seq1);
	if( length($seqL) > 0 ){
		my $endNew = $posD + len_leave($cigar_1, "MD") - 1;
		$seq_double_hash->{$chrD.":".$posD.":".$endNew.":".$cigar_1.":".$seqL} = "seqSame1";
	}
	($cigar_1,$cigar_2) = &sep_cigar($cigar_2, $#seq1, $#seq1);
	$posD = $endD - len_leave($cigar_2, "MD") + 1;
	my$seqR1 = substr($seqSeqD1, $n+1);
	my$seqR2 = substr($seqSeqD2, $n+1);
	if( $seqR1 eq $seqR2 and $seqR1 !~ /N/ ){
		$seq_double_hash->{$cigar_2.":".$seqR1} = "seqSame2";
	if( $seqR1 ne $seqR2 and length($seqR1)>1 ){
		$seq_double_hash = find_diff_seqSeq($chrD,$posD,$endD,$cigar_2, $seqR1, $seqR2, $seq_double_hash);
	}elsif(length($seqR1) == 1){
	}else{
		$seq_double_hash->{$cigar_2.":".$seqR1} = "seqSame2";
	}
	return($seq_double_hash);
}
sub sep_cigar{  #  分割cigar, M_pos 计入cigar1, 拼接时需注意
	my($matchInfo,$M_pos,$pos)=@_;	
	my($match_n_,$match_str_) = &sep_match_info($matchInfo);
	my @match_n = @$match_n_;
	my @match_str = @$match_str_;
	my ($seq_l,$seq_r,$seq_D,$seq_I,$seq_S,$seq_M)=(0,0,0,0,0,0);
	my ($cigar1, $cigar2, $store2) = ("","",0);
	foreach my $index(0..$#match_n){
		my $cigar = $match_n[$index].$match_str[$index];
		if( $match_str[$index] eq "D" ){
			if( $M_pos>=$seq_r+$pos and $M_pos<$match_n[$index]+$pos+$seq_r ){
				print "in delete $matchInfo !\n";
				return(0);
			}
			$seq_r += $match_n[$index];
			$seq_D += $match_n[$index];
		}elsif( $match_str[$index] eq "M" ){
			if( $M_pos>=$seq_r+$pos and $M_pos<$match_n[$index]+$pos+$seq_r ){
				$store2 = 1;
				$cigar = ($M_pos-$seq_r-$pos+1).$match_str[$index];
				$cigar1 .= $cigar;
				$cigar = ($match_n[$index]+$seq_r+$pos - $M_pos -1).$match_str[$index];
			}
			$seq_r += $match_n[$index]; 
			$seq_l += $match_n[$index];
			$seq_M += $match_n[$index];
		}elsif( $match_str[$index] eq "I" ){
			$seq_l += $match_n[$index];
			$seq_I += $match_n[$index];
		}elsif( $match_str[$index] eq "S" ){
			$seq_l += $match_n[$index];
			$seq_S += $match_n[$index];
		}
		
		if( $cigar =~ /^0/ ){ next; }
		if( $store2 ){	
			$cigar2 .= $cigar; 
		}else{
			$cigar1 .= $cigar;
		}
	}
	return($cigar1, $cigar2);
}
sub paste_cigar{  ### 2个cigar的粘贴
	my($cigar_a,$cigar_b)=@_;
	my ($match_n_a,$match_str_a) = &sep_match_info($cigar_a);
	my ($match_n_b,$match_str_b) = &sep_match_info($cigar_b);
	my ($re,$wrong);
	if( scalar(@$match_n_a) < 1 ){ return($cigar_b); }
	if( scalar(@$match_n_b) < 1 ){ return($cigar_a); }

	if( $match_str_a->[scalar(@$match_n_a)-1] eq $match_str_b->[0] ){
		foreach my$index(0..(scalar(@$match_n_a)-2)){ 	$re .= $match_n_a->[$index].$match_str_a->[$index]; }
		$re .= ($match_n_a->[scalar(@$match_n_a)-1] + $match_n_b->[0]).$match_str_b->[0];
		foreach my$index(1..(scalar(@$match_n_b)-1)){	$re .= $match_n_b->[$index].$match_str_b->[$index];	}
	}elsif( scalar(@$match_n_a) >1 and $match_str_b->[0] ne "S" and $match_str_a->[scalar(@$match_n_a)-1] eq "S" ){
		foreach my$index(0..(scalar(@$match_n_a)-2)){   $re .= $match_n_a->[$index].$match_str_a->[$index]; }
		$re .= ($match_n_a->[scalar(@$match_n_a)-1] + $match_n_b->[0])."S";
		foreach my$index(1..(scalar(@$match_n_b)-1)){   $re .= $match_n_b->[$index].$match_str_b->[$index]; }
		$wrong =1;
	}elsif( scalar(@$match_n_b) >1 and $match_str_b->[0] eq "S" and $match_str_a->[scalar(@$match_n_a)-1] ne "S" ){
		foreach my$index(0..(scalar(@$match_n_a)-2)){   $re .= $match_n_a->[$index].$match_str_a->[$index]; }
		$re .= ($match_n_a->[scalar(@$match_n_a)-1] + $match_n_b->[0])."S";
		foreach my$index(1..(scalar(@$match_n_b)-1)){   $re .= $match_n_b->[$index].$match_str_b->[$index]; }
		$wrong =1;
	}else{
		$re = $cigar_a.$cigar_b;
	}

	return($re);
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
sub explain_bam_flag{ # 解析 flags
	my ($flag)=@_;
	my $flag_bin=sprintf("%b", $flag);
	my @flag_bin = split //, $flag_bin;
	my $is_read1 = $flag_bin[-7];
	#my $is_read2 = @flag_bin>=8? $flag_bin[-8]: 0;
	#my $is_supplementary = @flag_bin>=12? $flag_bin[-12]: 0;
	#my $is_proper_pair = $flag_bin[-2];
	my $is_reverse = $flag_bin[-5];
	my $is_unmap = $flag_bin[-3];
	my $is_munmap = $flag_bin[-4];
	#my $dup = @flag_bin>=11? $flag_bin[-11]: 0;
	my $rnum = $is_read1==1? 1: 2;
	#return($rnum, $is_proper_pair, $is_reverse, $is_unmap, $is_munmap, $is_supplementary);
	return($rnum,$is_reverse,$is_unmap);
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

