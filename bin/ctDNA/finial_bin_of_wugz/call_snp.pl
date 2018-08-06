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
my ($sampleid,$bamfile,$outdir,$regionfile,$regionprimerfile,$outfile);
my ($read1only,$endborderN,$print_type);
my ($testchr,$testpos);
GetOptions(
				"help|?" =>\&USAGE,
				"s:s"=>\$sampleid,
				"i:s"=>\$bamfile,
				"o:s"=>\$outfile,
				"e:s"=>\$endborderN,
				"r1:s"=>\$read1only,
				"d:s"=>\$outdir,
				"r:s"=>\$regionfile,
				"rp:s"=>\$regionprimerfile,
				"p:s"=>\$print_type,
				) or &USAGE;
&USAGE unless defined($bamfile and $outfile);

if(not defined $print_type){	$print_type=0;	}
if(not defined $endborderN){	$endborderN=0;	}
if(not defined $read1only){		$read1only=0; }
$sampleid||="Sample";
$testchr||="chr11";
$testpos||=5247905;
$outdir||="./";
`mkdir $outdir` unless (-d $outdir);
$outdir=&AbsolutePath("dir",$outdir);
open(I, "samtools view $bamfile|") or die $!;
open(OUT,">$outfile");
my %hash;	# {$id}{$chr}{"READ$rnum"}
my %count;	# {$chr}{$pos}{$primer}{ substr($seq2,$i,1) }
my %count_; my %count_z; 
my %id_primer;	# $id_primer{$id} = $primer
my %id_fangxiang;
my %ref_del;	# {$chr}{$pos}{$primer}{ $delete }
my %ref_del_; my %ref_del_z; 
my %ref_insert;	# {$chr}{$pos}{$primer}{ $insert }
my %ref_insert_; my %ref_insert_z; 
my %ref_hash;	# {$chr}{$pos} 存储存在突变位置的
my %ref_hash_del;   # {$id}{$chr}{$pos}{$pos+$match_n[$index]-1}
my %rnum_hash;	# 1=>2 ; 2=>1 ;
my %paired_Snp; # {$id}{"READ".$rnum} = [$chr,$pos,$pos,$primer,$insert]
				# 此处暂未优化，可进一步优化，思路类似 indel
my %paired_Del;	# {$id}{"READ".$rnum} = [$chr,$pos,$pos+$len-1,$primer,$delete]
my %paired_In;	# {$id}{"READ".$rnum} = [$chr,$pos-1,$pos,$primer,$insert]
my $RefGenome = "/data/bioit/biodata/duyp/bin/hg19/hg19.fasta";

&load_paired_rnum();
&print_file_header($sampleid);
&SHOW_TIME("Reading");
while(<I>){
	chomp; my$lines = $_;
	&check_and_storage_matchInfo($lines);
}
close I;
&SHOW_TIME("Counting");
&check_and_delete_overlap();
&SHOW_TIME("printing:");
my (%snp_count,%snp_count_2,%snp_count_z,%snp_count_f);
if(defined $regionfile){
	my %position_print;
	open(R,"$regionfile")or die $!;
	while(<R>){
		my($r_chr,$r_s,$r_e,@line) = split(/\t/,$_);
		foreach my $pos( $r_s..$r_e ){
			$position_print{$r_chr}{$pos} = 1;
		}
	}
	close R;

	foreach my $chr(sort keys%position_print){
		my$position_print_2 = $position_print{$chr};
		foreach my $pos( sort{ $b <=> $a } keys %$position_print_2){
			my $hash2chr = $count{$chr};
			if(exists $hash2chr->{$pos}){
				&print_base($chr,$pos,$hash2chr->{$pos}, $count_{$chr}{$pos},$count_z{$chr}{$pos},"snp");
			}	
			if( exists $ref_insert{$chr}{$pos} ){
				&print_base($chr,$pos,$ref_insert{$chr}{$pos},$ref_insert_{$chr}{$pos},$ref_insert_z{$chr}{$pos},"insert");
			}elsif( exists $ref_del{$chr}{$pos} ){
				&print_base($chr,$pos,$ref_del{$chr}{$pos},$ref_del_{$chr}{$pos},$ref_del_z{$chr}{$pos},"delete");
			}
		}
	}
}elsif(defined $regionprimerfile){
	my %position_print;
	open(R,"$regionprimerfile")or die $!;
	while(<R>){
		my($r_chr,$r_s,$r_e,$special_primer) = split(/\t/,$_);
		my @a_primer = split /,/,$special_primer;
		foreach my $pos( $r_s..$r_e ){
			@{$position_print{$r_chr}{$pos}} =@a_primer;
		}
	}
	close R;

	foreach my $chr(sort keys%position_print){
		my$position_print_2 = $position_print{$chr};
		foreach my $pos( sort{ $b <=> $a } keys %$position_print_2){
			my $hash2chr = $count{$chr}; #&together($hash2chr->{$pos},$position_print{$chr}{$pos})
			my @primer = @{$position_print{$chr}{$pos}};
			if(exists $hash2chr->{$pos} ){ 
				&print_base($chr,$pos,$hash2chr->{$pos}, $count_{$chr}{$pos},$count_z{$chr}{$pos},"snp",\@primer);
			}
			if( exists $ref_insert{$chr}{$pos} ){
				&print_base($chr,$pos,$ref_insert{$chr}{$pos},$ref_insert_{$chr}{$pos},$ref_insert_z{$chr}{$pos},"insert",\@primer);
			}elsif( exists $ref_del{$chr}{$pos} ){
				&print_base($chr,$pos,$ref_del{$chr}{$pos},$ref_del_{$chr}{$pos},$ref_del_z{$chr}{$pos},"delete",\@primer);
			}
		}
	}	
}else{
	foreach my$chr(sort keys %count){
		my $hash2chr = $count{$chr};
		foreach my $pos (sort{ $b <=> $a } keys %$hash2chr){
			&print_base($chr,$pos,$hash2chr->{$pos},$count_{$chr}{$pos},$count_z{$chr}{$pos},"snp");
			if( exists $ref_insert{$chr}{$pos} ){
			#	print "\n$chr,$pos,$ref_snp\tinsert???\n";
				&print_base($chr,$pos,$ref_insert{$chr}{$pos},$ref_insert_{$chr}{$pos},$ref_insert_z{$chr}{$pos},"insert");
			}elsif( exists $ref_del{$chr}{$pos} ){
				&print_base($chr,$pos,$ref_del{$chr}{$pos},$ref_del_{$chr}{$pos},$ref_del_z{$chr}{$pos},"delete");
			}
		}
	}
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

&SHOW_TIME("Ending");
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub load_paired_rnum{
	$rnum_hash{"1"}="2";
	$rnum_hash{"2"}="1";
}
sub together{
	my($a,$b) = @_;
	my %h;
	foreach(@$b){ $h{$_}=1; }
	foreach(keys %$a){
		if( exists $h{$_} ){ return 1; }
	}
	return 0;
}

sub print_base{
	my( $chr ,$pos ,$hash3pos,$hash3pos_,$hash3pos_z,$snp,$special_primer,$print_all_DP ) = @_;
	my (@primer_,$PR_re,$ADPR_re);
	if(not defined $print_all_DP){
		($PR_re,$ADPR_re) = &print_base( $chr ,$pos ,$hash3pos,$hash3pos_,$hash3pos_z,$snp,$special_primer,"1" );
		if(defined $special_primer){ 
			@primer_ = sort @$special_primer; 
		}else{
			@primer_ = sort keys %$hash3pos;
		}
	}else{
		@primer_ = sort keys %$hash3pos;
	}
	
	my $printout;
	if( $snp ne "delete" ){
		$printout = "$chr\t$pos\t.\t";
	}else{
		$printout = "$chr\t".($pos-1)."\t.\t";
	}
	my ($count_all,$DP2,$DP1,$PR,$ADPR) = (0,0,0,"PR=","ADPR=") ;
	my (%base_N_hash,%f_N_hash,$reference,$alt,$alt_min_l);

	foreach my $primer ( sort @primer_ ){
		my $hash4primer = $hash3pos->{$primer};
		foreach my $base_N (sort keys %$hash4primer){ 
			if( $base_N eq "N" ){next;}
			( exists $base_N_hash{$base_N} )?($base_N_hash{$base_N} = $base_N_hash{$base_N} + $hash4primer->{$base_N}):
				($base_N_hash{$base_N} = $hash4primer->{$base_N});
		}
	}

	if (%base_N_hash){}else{ return; }
	if(exists $ref_hash{$chr}{$pos} and $snp eq "snp"){
		$reference = $ref_hash{$chr}{$pos};	
		foreach my $base_N(sort keys%base_N_hash){
			if( $base_N ne $ref_hash{$chr}{$pos}){
				$alt .= "$base_N,";
			}
		}
		$alt .= "<*>";
#		if(defined $alt){ 
#			$alt .= "<*>"; 
#		}else{ 
#			print "$chr,$pos,$snp  is in trouble!\t$printout\n"; 
#			return; 
#		}
	}elsif( $snp eq "snp" ){
		$reference = join '',keys%base_N_hash;
		$alt .= "<*>";
	}elsif( $snp eq "insert" or $snp eq "delete" ){
		foreach my $base_N(sort keys%base_N_hash){
			if( $base_N eq "N" ){next;}
			if(defined $reference){ 
				my($reference_,$indel_)=&check_indel($chr,$pos+1,$base_N,$snp);
				if( length($reference_) >= length($reference) ){
					$reference = $reference_;
				}
			}else{
				my($reference_,$indel_)=&check_indel($chr,$pos+1,$base_N,$snp);
				$reference = $reference_;
				#$alt .= $indel_.",";
			}
		}
		foreach my $base_N(sort keys%base_N_hash){
			if( $base_N eq "N" ){next;}
			if( $snp eq "insert" ){
				my $indel_ = $reference.$base_N;
				$alt .= $indel_.",";
			}else{
				if( not defined $alt_min_l ){	$alt_min_l = length($base_N);	}
				if( defined $alt_min_l and $alt_min_l > length($base_N) ){	$alt_min_l = length($base_N);	}
				my $indel_ = $reference;
				$indel_ =~ s/$base_N//;
				$alt .= $indel_.",";
			}
		}
		chop($alt);
	}
	if(not defined $reference or not defined $alt){	print keys%base_N_hash,"##$reference\t$alt\tnot defined ?\n";	}
	$printout .= "$reference\t$alt\t0\t.\t";
	if( $snp ne "snp" ){ $printout .= "INDEL;"; }

	($count_all,$DP2,$DP1,$PR,$ADPR) = (0,0,0,"PR=","ADPR=") ;
	%base_N_hash=();%f_N_hash=();	
	foreach my $primer ( @primer_ ){
		$PR .= $primer.",";
		my $count_primer = 0 ;
		my $hash4primer = $hash3pos->{$primer};
		if( $snp eq "snp" ){
			my $base_N = $reference ;
			if( exists $hash4primer->{$base_N} ){
			}else{
				$hash3pos_z->{$primer}->{$base_N} = 0;
				$hash4primer->{$base_N} = 0;
				$hash3pos_->{$primer}->{$base_N} = 0;
			}
			( exists $f_N_hash{$base_N} )?
				($f_N_hash{$base_N} += $hash3pos_z->{$primer}->{$base_N}):
				($f_N_hash{$base_N}  = $hash3pos_z->{$primer}->{$base_N});
			( exists $base_N_hash{$base_N} )?
				($base_N_hash{$base_N} = $base_N_hash{$base_N} + $hash4primer->{$base_N}):
				($base_N_hash{$base_N} = $hash4primer->{$base_N});
#			$ADPR .= $hash4primer->{$base_N}.",";
			$ADPR .= $base_N."-".$hash4primer->{$base_N}.",";
			$count_all += $hash4primer->{$base_N};
			$DP1 += $hash4primer->{$base_N};
			$DP2 += $hash3pos_->{$primer}->{$base_N};
#			$ADR += $hash3pos_z->{$primer}->{$base_N} ;
		}
		foreach my $base_N (sort keys %$hash4primer){
			if( $snp eq "snp" and $base_N eq $reference ){next;}
			if( $base_N eq "N" ){next;}  ##################################
			# 此处若有N insert，则导致 p 偏低
			# 此处若有N snp，则导致 p 偏低
			( exists $f_N_hash{$base_N} )?
				($f_N_hash{$base_N} += $hash3pos_z->{$primer}->{$base_N}):
				($f_N_hash{$base_N}  = $hash3pos_z->{$primer}->{$base_N});
			( exists $base_N_hash{$base_N} )?
				($base_N_hash{$base_N} = $base_N_hash{$base_N} + $hash4primer->{$base_N}):
				($base_N_hash{$base_N} = $hash4primer->{$base_N});
#			$ADPR .= $hash4primer->{$base_N}.",";
			$ADPR .= $base_N."-".$hash4primer->{$base_N}.","; 
			$count_all += $hash4primer->{$base_N};
			$DP1 += $hash4primer->{$base_N};
			$DP2 += $hash3pos_->{$primer}->{$base_N};
#			$ADR += $hash3pos_z->{$primer}->{$base_N} ;
		}
		chop($ADPR);$ADPR .= ":";
	}
	if($print_all_DP){
		return($PR,$ADPR);	
	}else{
		$PR = $PR_re;
		$ADPR = $ADPR_re;
	}
	chop($PR);
	$printout .= "DP=$DP2;$PR;";
	chop($ADPR);
	$printout .= "$ADPR;";
	$printout .= "AD=";  
	my ($printAD,$printADF,$printADR,$printDP,$ADF,$ADR);
	my ($ad,$adf,$adr,$dp) =(0,0,0,0);
	if( $snp eq "snp" ){
		if(exists $base_N_hash{$reference}){
			$printAD = $base_N_hash{$reference};
#			$printDP = $DP2;    changed this DP to AD
			$printDP = $count_all;
			$printADF = $base_N_hash{$reference} - $f_N_hash{$reference};
			$printADR = $f_N_hash{$reference};
		}
		foreach my$key(sort keys%base_N_hash){ 
			$snp_count{$chr}{$pos} += $base_N_hash{$key}; 
			$snp_count_z{$chr}{$pos} += $base_N_hash{$key} - $f_N_hash{$key};
			$snp_count_f{$chr}{$pos} += $f_N_hash{$key};
			if( $key eq $reference ){next;}
			$printAD .=  ",".$base_N_hash{$key};
			$ADF = $base_N_hash{$key} - $f_N_hash{$key};
			$printADF .=  ",".$ADF;
			$ADR = $f_N_hash{$key};
			$printADR .=  ",".$ADR;
		}
#		$snp_count_2{$chr}{$pos} = $DP2;    changed this DP to AD
		$snp_count_2{$chr}{$pos} = $count_all;
		$printAD .=  ",0";
		$printADF .=  ",0";
		$printADR .=  ",0"
	}elsif( $snp eq "delete" ){
		foreach my$key(sort keys%base_N_hash){
			$printAD .=  ",".$base_N_hash{$key};
			$ADF = $base_N_hash{$key} - $f_N_hash{$key};
			$printADF .=  ",".$ADF;
			$ADR = $f_N_hash{$key};
			$printADR .=  ",".$ADR;
		}
		if( not defined $snp_count_2{$chr}{$pos} ){ 
			$snp_count_2{$chr}{$pos} = 0; 
			$snp_count{$chr}{$pos} = 0;
			$snp_count_z{$chr}{$pos} = 0;
			$snp_count_f{$chr}{$pos} = 0;
		}
		if( not defined $snp_count_2{$chr}{$pos+$alt_min_l-1} ){
			$snp_count_2{$chr}{$pos+$alt_min_l-1} = 0;
			$snp_count{$chr}{$pos+$alt_min_l-1} = 0;
			$snp_count_z{$chr}{$pos+$alt_min_l-1} = 0;
			$snp_count_z{$chr}{$pos+$alt_min_l-1} = 0;
		}
#		$dp = $DP2+$snp_count_2{$chr}{$pos};   changed this DP to AD
		$dp = $count_all + $snp_count_2{$chr}{$pos};
		$ad = $snp_count{$chr}{$pos};
		$adf = $snp_count_z{$chr}{$pos};
		$adr = $snp_count_f{$chr}{$pos};
		if( $pos < $pos+$alt_min_l-1 and $snp_count{$chr}{$pos} < $snp_count{$chr}{$pos+$alt_min_l-1} ){
			$dp = $count_all + $snp_count_2{$chr}{$pos+$alt_min_l-1};
			$ad = $snp_count{$chr}{$pos+$alt_min_l-1};
			$adf = $snp_count_z{$chr}{$pos+$alt_min_l-1};
			$adr = $snp_count_z{$chr}{$pos+$alt_min_l-1};
		}
		$printAD = $ad.$printAD;
		$printDP = $dp;
		$printADF = $adf.$printADF;
		$printADR = $adr.$printADR; 
	}elsif( $snp eq "insert" ){
		foreach my$key(sort keys%base_N_hash){
			$printAD .=  ",".$base_N_hash{$key};
			$ad += $base_N_hash{$key};
			$ADF = $base_N_hash{$key} - $f_N_hash{$key};
			$printADF .=  ",".$ADF;
			$adf += $base_N_hash{$key} - $f_N_hash{$key};
			$ADR = $f_N_hash{$key};
			$printADR .=  ",".$ADR;
			$adr += $f_N_hash{$key};
		}
		if( not defined $snp_count_2{$chr}{$pos} ){ 
			$snp_count_2{$chr}{$pos} = 0; 
			$snp_count{$chr}{$pos} = 0;
			$snp_count_z{$chr}{$pos} = 0;
			$snp_count_f{$chr}{$pos} = 0;
		}
		if( not defined $snp_count_2{$chr}{$pos+1} ){ 
			$snp_count_2{$chr}{$pos+1} = 0; 
			$snp_count{$chr}{$pos+1} = 0;
			$snp_count_z{$chr}{$pos+1} = 0;
			$snp_count_f{$chr}{$pos+1} = 0;
		}
		$dp = $snp_count_2{$chr}{$pos};
		$ad = $snp_count{$chr}{$pos} - $ad;
		$adf = $snp_count_z{$chr}{$pos} - $adf;
		$adr = $snp_count_f{$chr}{$pos} - $adr;
		if( $snp_count{$chr}{$pos+1} > $snp_count{$chr}{$pos} ){
			$dp = $snp_count_2{$chr}{$pos+1};
			$ad = $snp_count{$chr}{$pos+1} - $ad;
			$adf = $snp_count_z{$chr}{$pos+1} - $adf;
			$adr = $snp_count_f{$chr}{$pos+1} - $adr;		}
		$printAD = $ad.$printAD;
		$printDP = $dp;
		$printADF = $adf.$printADF;
		$printADR = $adr.$printADR;
		if( $ad < 0 ){
			#print "$chr,$pos,$snp  is in trouble!\t$printout\n";	
 			return;
		}	
	}
	$printout .= "$printAD\tDP:ADF:ADR:AD\t$printDP:$printADF:$printADR:$printAD";
	if( $printDP <= 0 ){
		#print "$chr,$pos,$snp  is in trouble!\t$printout\n";
		return;
	}
	print OUT $printout,"\n";
}

sub SHOW_TIME {
	#显示当时时间函数，参数内容是时间前得提示信息，为字符串
	my ($str)=@_;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	my $temp=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
	print "$str:\t[".$temp."]\n";
}
sub explain_bam_flag{
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
sub sep_match_info{
	my($match_info)=@_;
	my @ucigar = split //, $match_info;
	my (@match_n,@match_str);
	my $cigar="";
	foreach my $i(0..$#ucigar){
		if($ucigar[$i] eq "M" || $ucigar[$i] eq "I" || $ucigar[$i] eq "D" || $ucigar[$i] eq "S"){
			push @match_str, $ucigar[$i];
			push @match_n, $cigar;
			$cigar = "";
		}else{
			$cigar .= $ucigar[$i];
		}
	}
	return(\@match_n,\@match_str)
}
sub sep_MD_info{
	my($ref_or_not,$chr,$pos,$id)=@_;
	my @md = split //, $ref_or_not;
	my $n_pos = 0; my $md_n= 0;
	foreach my $i(0..$#md){
		if($md[$i] =~ /\d/){
			$md_n .= $md[$i];
		}elsif($md[$i] =~ /([ATGCN])/){
			my $base_ref = $1;
			#print "$ref_or_not\t$chr\t$pos\t$n_pos,$md_n\n";
			$n_pos += $md_n ; # 加上每次正确匹配长度，即此点reads上前端长度
			$ref_hash{$chr}{$n_pos + $pos} = $base_ref;
			$n_pos++; # 记录一个基因组位点，相对位置加1，即前端长度加1
			$md_n= 0;
		}elsif($md[$i] =~ /\^/){
			my $del_base_ref = "";
			foreach my$j(1..$#md){
				if($md[$i+$j] =~ /([ATGCN])/){
					$del_base_ref .= $1;
				}else{ last; }
			}
			$ref_hash_del{$id}{$chr}{$n_pos+$md_n+$pos}{$n_pos+$md_n+$pos+length($del_base_ref)-1} = $del_base_ref;
			#$i = $i+$j;
		}
	}
}
sub seq_l_endborderN{
	my($index,$seq_l_start,$seq_l_end,$pos,$len,$read_pos,$lenseq,$maxindex,$match_str) = @_;
	my $lenMI = $len;
    if( $match_str eq "D" ){ $lenMI = 0; }
	my($posM,$lenM,$read_posM,$ID_next) = ($pos,$len,$read_pos,0);
	my $len_leave = $seq_l_end - $seq_l_start;
	if( $seq_l_start > 1 + $read_pos ){
		$ID_next = 1;
		$posM = $pos + $seq_l_start-1 - $read_pos;
		$lenM = $len - $seq_l_start+1 + $read_pos;
		$read_posM = $seq_l_start - 1 ;
		if( $seq_l_end - ($read_pos + $lenMI ) < 0 ){
			$lenM = $lenM - ($read_pos + $lenMI + 1 - $seq_l_end);
		}
	}elsif( $seq_l_end - ($read_pos + $lenMI ) < 0 ){
		$ID_next = 1;
		$lenM = $len - ($read_pos + $lenMI - $seq_l_end);
	} #print "$posM,$lenM,$read_posM,$ID_next\t$seq_l_start\t$read_pos+1\t$seq_l_end - ($read_pos + $lenMI)\n";
	return($posM,$lenM,$read_posM,$ID_next);
}
sub sep_match_info_count_ref{
	my($chr,$pos,$id,$rnum,$primer,$seq,$seq_ref,$match_n_,$match_str_,$match_info)=@_;
	my@match_n = @$match_n_; my @match_str = @$match_str_;
	my ($start_pos,$read_pos,$read_ref_pos) = ($pos, 0, 0);
	my ($seq_l_start, $seq_l_end) = (1, length($seq));
	if( $id_fangxiang{$id}{$rnum} ){ 
		$seq_l_start = $seq_l_start + $endborderN; } # 去末端
	else{	
		$seq_l_end = $seq_l_end - $endborderN;	}	# 去末端
	if( $rnum == 2 and $id_fangxiang{$id}{$rnum} == 1){ $seq_l_end = $seq_l_end - $endborderN; } # 去read2始端
	if( $rnum == 2 and $id_fangxiang{$id}{$rnum} == 0){ $seq_l_start = $seq_l_start + $endborderN; } # 去read2始端
	if( $match_str[0] eq "S" ){	$seq_l_start = $seq_l_start + $match_n[0];	}				# 去S
	if( $match_str[$#match_n] eq "S" ){ $seq_l_end = $seq_l_end - $match_n[$#match_n];  }
	if($print_type == 2){
		print "$id,read$rnum\t$seq_l_start, $seq_l_end\t$id_fangxiang{$id}{$rnum}\t#@match_n#@match_str\n";
	}
	foreach my $index(0..$#match_n){
		if($match_str[$index] eq "S" && $index!=0 && $index!=$#match_n){# 若存在 S 在中间则报错
			print "Wrong: $match_info 'S' is in middle!\n";die; }
		my $len = $match_n[$index];
		if($match_str[$index] eq "M"){
			my($posM,$lenM,$read_posM,$ID_next)=
				&seq_l_endborderN($index,$seq_l_start,$seq_l_end,$pos,$len,$read_pos,length($seq),$#match_n,$match_str[$index]);
			if($print_type == 2){
				print substr($seq,$read_posM,$lenM),"\t";
				if( $seq_ref ne "noref" ){	print substr($seq_ref,$posM-$start_pos,$lenM),"\t";	}
			print "what:$id,$rnum:$pos,$len,$read_pos\t$posM,$lenM,$read_posM,$ID_next\t$match_n[$index],$match_str[$index]\n";}
			if( $lenM <= 0 ){	
				$pos += $len;  $read_pos += $len;  $read_ref_pos += $len;
				next;	}
			if( $seq_ref eq "noref" ){  
				&countAdd($chr,$posM,$lenM,$id,$rnum,$primer, substr($seq,$read_posM,$lenM) ); #+1
			}else{
				&countAdd_ref($chr,$posM,$lenM,$id,$rnum,$primer, substr($seq,$read_posM,$lenM),substr($seq_ref,$posM-$start_pos,$lenM)); #+1
			}
			push @{$hash{$id}{$chr}{"READ".$rnum}},[$posM, $posM+$lenM-1 ,substr($seq,$read_posM,$lenM)];
			$pos += $len;  $read_pos += $len;  $read_ref_pos += $len;
		}elsif($match_str[$index] eq "D"){
			my($posM,$lenM,$read_posM,$ID_next)=
				&seq_l_endborderN($index,$seq_l_start,$seq_l_end,$pos,$len,$read_pos,length($seq),$#match_n,$match_str[$index]);
			if( $ID_next ){	
				$pos = $pos+$len; $read_ref_pos += $len;
				if($print_type == 2){
			print "what:$id,$rnum:$pos,$len,$read_pos\t$posM,$lenM,$read_posM,$ID_next\t$match_n[$index],$match_str[$index]\n";	}
				next;	}
			my $delete = "no";
			if( $seq_ref eq "noref" ){ 
				$delete = $ref_hash_del{$id}{$chr}{$pos}{$pos+$match_n[$index]-1};		
			}else{
				$delete =  substr($seq_ref,$read_ref_pos,$len);
			}	
#			if( $id eq "1322-2712_8" ){ print "$delete\t$id,$chr,$pos\t$pos+$match_n[$index]-1\t$ref_hash_del{$id}{$chr}{$pos}{$pos+$match_n[$index]-1}\n"; }

			if( exists $paired_Del{$id}{"READ".$rnum_hash{$rnum}}{$pos} ){
				my ($chr_,$start_,$end_,$primer_,$seq_) = @{$paired_Del{$id}{"READ".$rnum_hash{$rnum}}{$pos}};
				if( $chr eq $chr_ and $pos eq $start_ and $pos+$len-1 eq $end_ and $seq_ eq $delete ){ 
					$ref_del_z{$chr}{$pos}{$primer}{ $delete } += $id_fangxiang{$id}{$rnum};
				}else{ 
					print "$id,READ$rnum,$chr,$pos\tat this position, delete of read1 & read2 are not same!\n";
					print "$chr_,$start_,$end_,$primer_,$seq_\tREAD$rnum_hash{$rnum}\n";
					print "$chr,$pos,$pos+$len-1,$delete\tREAD$rnum\n";
					$ref_del{$chr_}{$start_}{$primer_}{$seq_}--;
					if( $ref_del{$chr_}{$start_}{$primer_}{$seq_} == 0 ){
						delete $ref_del{$chr_}{$start_}{$primer_}{$seq_};}
					$ref_del_z{$chr_}{$start_}{$primer_}{$seq_} -= $id_fangxiang{$id}{$rnum}; 
				}
			}else{
				$ref_del{$chr}{$pos}{$primer}{ $delete }++; #存del
				$ref_del_z{$chr}{$pos}{$primer}{ $delete } += $id_fangxiang{$id}{$rnum};
				$paired_Del{$id}{"READ".$rnum}{$pos} = [$chr,$pos,$pos+$len-1,$primer,$delete];
				#print "$id,READ$rnum,$chr,$pos\tat this position, delete of read1 !\n";
			}
			$ref_del_{$chr}{$pos}{$primer}{ $delete }++;
			$pos = $pos+$len; $read_ref_pos += $len;
		}elsif($match_str[$index] eq "S" ){
			$read_pos += $len;
		}elsif( $match_str[$index] eq "I" ){
			my($posM,$lenM,$read_posM,$ID_next)=
				&seq_l_endborderN($index,$seq_l_start,$seq_l_end,$pos,$len,$read_pos,length($seq),$#match_str,$match_str[$index]);
			if( $ID_next ){ 
				$read_pos += $len;
				if($print_type == 2){
			print "what:$id,$rnum:$pos,$len,$read_pos\t$posM,$lenM,$read_posM,$ID_next\t$match_n[$index],$match_str[$index]\n";	}
				next;   }
			my $insert = substr($seq,$read_pos,$len);
			if( exists $paired_In{$id}{"READ".$rnum_hash{$rnum}}{$pos-1} ){
				my ($chr_,$start_,$end_,$primer_,$seq_) = @{$paired_In{$id}{"READ".$rnum_hash{$rnum}}{$pos-1}};
				if( $chr eq $chr_ and $pos-1 eq $start_ and $pos eq $end_  and $seq_ eq $insert ){
					$ref_insert_z{$chr}{$pos-1}{$primer}{ $insert } += $id_fangxiang{$id}{$rnum}; 
				}else{
					print "$id,$chr,$pos-1\tat this position, insert of read1 & read2 are not same!\n";
					$ref_insert{$chr_}{$start_}{$primer_}{$seq_}--;
					if( $ref_insert{$chr_}{$start_}{$primer_}{$seq_} == 0 ){
						delete $ref_insert{$chr_}{$start_}{$primer_}{$seq_}; }
					$ref_insert_z{$chr_}{$start_}{$primer_}{$seq_} -= $id_fangxiang{$id}{$rnum};
				}
			}else{
				$ref_insert{$chr}{$pos-1}{$primer}{ $insert }++;
				$ref_insert_z{$chr}{$pos-1}{$primer}{ $insert } += $id_fangxiang{$id}{$rnum};
				$paired_In{$id}{"READ".$rnum}{$pos-1} = [$chr,$pos-1,$pos,$primer,$insert];
			}
			$ref_insert_{$chr}{$pos-1}{$primer}{ $insert }++;
			$read_pos += $len;
		}
	}
}

sub sep_and_count_array{
	my($match_info,$chr,$pos,$id,$rnum,$primer,$seq,$seq_ref) = @_;
	my ($match_n,$match_str) = &sep_match_info($match_info);# 解析 match_info
	&sep_match_info_count_ref($chr,$pos,$id,$rnum,$primer,$seq,$seq_ref,\@$match_n,\@$match_str,$match_info);
	
}

sub check_and_storage_matchInfo{
	my ($lines) = @_ ;
	my ($primer,$ref_or_not );
	my ($id,$flag,$chr,$pos,$mapQ,$matchInfo,$li7,$li8,$li9,$seq,@line) = split(/\t/,$lines);
	my($rnum, $is_reverse ,$is_unmap)=&explain_bam_flag($flag);
	if( $read1only == 1 and $rnum == 2){ return; }

	if( $lines =~ m/XP:Z:(.+?)\t/){ 
		$primer = $1; 
		$id_primer{$id} = $primer;
		$id_fangxiang{$id}{$rnum} = $is_reverse;
	}else{ 
		print "no primer?$lines\n";die; 
	}

	if( $lines =~ m/MD:Z:(.+?)\t/){ 
		$ref_or_not = $1;
		&sep_MD_info($ref_or_not,$chr,$pos,$id); # 解析 ref_or_not，存储突变位点的参考基因组位置 
		if($is_unmap == 1){ #   print "$id,$flag,$chr,$pos,unmap\n";
		}else{
			&sep_and_count_array($matchInfo,$chr,$pos,$id,$rnum,$primer,$seq,"noref");
		}
	}else{ 	#print "no MD?\t$lines\n";
		my $end = $pos+100;
		my $seq_ref = `samtools faidx /data/bioit/biodata/duyp/bin/hg19/hg19.fasta $chr:$pos-$end`;
		$seq_ref  =~ tr/ATGCatgc/ATGCATGC/; # 无MD，从基因组获取 seq_ref
		my( $ref_id, $seq_ref_1, $seq_ref_2 ) = split( /\n/,$seq_ref);
		if($is_unmap == 1){ #   print "$id,$flag,$chr,$pos,unmap\n";
		}else{
			&sep_and_count_array($matchInfo,$chr,$pos,$id,$rnum,$primer,$seq,$seq_ref_1.$seq_ref_2);
		}
	}
}

sub check_and_delete_overlap{
	foreach my$id(keys %hash){
		my $primer = $id_primer{$id};
		foreach my $chr(keys %{$hash{$id}}){
			if( exists $hash{$id}{$chr}{"READ1"} and exists $hash{$id}{$chr}{"READ2"} ){ 
				my @read1 = @{$hash{$id}{$chr}{"READ1"}};
				my @read2 = @{$hash{$id}{$chr}{"READ2"}};
				for(my $i=0; $i<@read1; $i++){	
					my ($start1,$end1,$seq1) = @{$read1[$i]};
					for (my $j=0; $j<@read2; $j++){
						my ($start2,$end2,$seq2) = @{$read2[$j]};
						if( $start1 <= $start2 ){
							if( $start2<=$end1 and $end1<=$end2 ){ # start1 start2 end1 end2
								&countDel( $chr ,$start2 ,$end1-$start2+1 ,$id,$primer,
								substr($seq1,$start2-$start1,$end1-$start2+1),
								substr($seq2,0,$end1-$start2+1)  );
							#&check_info_print("checkregion",$chr1,$start2,$end1,$k1,"del");
							}elsif( $start2<=$end1 and $end1>$end2 ){ # start1 start2 end2 end1
								&countDel( $chr ,$start2 ,$end2-$start2+1 ,$id,$primer,
								substr($seq1,$start2-$start1,$end2-$start2+1), 
								substr($seq2,0,$end2-$start2+1) );
							#&check_info_print("checkregion",$chr1,$start2,$end2,$k1,"del");
							}
						}else{
							if( $start1<=$end2 and $end2<=$end1 ){ # start2 start1 end2 end1
								&countDel( $chr ,$start1 ,$end2-$start1+1 ,$id,$primer,
								substr($seq2,$start1-$start2,$end2-$start1+1),
								substr($seq1,0,$end2-$start1+1)  );
							#&check_info_print("checkregion",$chr1,$start1,$end2,$k1,"del");
							}elsif( $start1<=$end2 and $end2>$end1 ){ # start2 start1 end1 end2
								&countDel( $chr ,$start1 ,$end1-$start1+1 ,$id,$primer,
								substr($seq2,$start1-$start2,$end1-$start1+1),
								substr($seq1,0,$end1-$start1+1)  );
							#&check_info_print("checkregion",$chr1,$start1,$end1,$k1,"del");
							}
						}
					}
				}
			}
		}
	}
}


sub AbsolutePath{		#获取指定目录或文件的决定路径
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
sub check_info_print{
	my($checktype,$chr,$pos,$len,$type,$j)=@_;
	if($checktype eq "checkpos"){
		if( $chr eq $testchr and $j == $testpos ){
			print "$checktype,$chr,$pos,$len,$type,$j\n";
		}
	}elsif($checktype eq "checkregion"){
		if( $chr eq $testchr and $testpos >=$pos and $testpos <= $len ){
			print "$checktype,$chr,$pos,$len,$type,$j\n";
		}
	}
}
sub countAdd{
	my($chr,$pos,$len,$id,$rnum,$primer,$seq)=@_;
	foreach my$i(0..(length($seq)-1)){
		#&check_info_print("checkpos",$chr,$pos,$len,$type,$j);
		$count{$chr}{$pos+$i}{$primer}{ substr($seq,$i,1) }++;
		$count_{$chr}{$pos+$i}{$primer}{ substr($seq,$i,1) }++;
		$count_z{$chr}{$pos+$i}{$primer}{ substr($seq,$i,1) } += $id_fangxiang{$id}{$rnum};
	}
}
sub countAdd_ref{
	my($chr,$pos,$len,$id,$rnum,$primer,$seq1,$seq2)=@_;
	if($seq1 eq $seq2){
		foreach my$i(0..(length($seq1)-1) ){
			$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }++;
			$count_{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }++;
			$count_z{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } += $id_fangxiang{$id}{$rnum};
		}
	}else{
		foreach my$i(0..(length($seq1)-1) ){
			$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }++;
			$count_{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }++;
			$count_z{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } += $id_fangxiang{$id}{$rnum};
			#print "read#$seq1\tref#$seq2\t$i\n";
			if( substr($seq1,$i,1) ne substr($seq2,$i,1) ){
				$ref_hash{$chr}{$pos+$i} = substr($seq2,$i,1);
			}
		}
	}
}
sub countDel{
	my($chr,$pos,$len,$id,$primer,$seq1,$seq2) = @_;
	if($seq1 eq $seq2){
		foreach my$i(0..(length($seq1)-1) ){
			if(not exists$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } ){
				print "$chr,$pos,$len,$i,$primer:",substr($seq1,$i,1),"\tnotadd\n"; next;}
			$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }--;
#			$count_z{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } -= $id_fangxiang{$id};
			if( $count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } == 0 ){ 
				delete $count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }; }
		} #print "del sucess\n";
	}else{#print "del false \t$seq1\t$seq2\n";
		foreach my$i(0..(length($seq1)-1) ){
			if( substr($seq1,$i,1) eq substr($seq2,$i,1) ){
				if(not exists$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }){
					print "$chr,$pos,$len,$i,$primer:",substr($seq1,$i,1),"\tnotadd\n"; next;}
				$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }--;
#				$count_z{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } -= $id_fangxiang{$id}{$rnum};
				if( $count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } == 0 ){
					delete $count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }; }
			}else{
				if(not exists$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }){
					print "$chr,$pos,$len,$i,$primer:",substr($seq1,$i,1),"\tnotadd\n"; next;}
				$count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }--;
				$count_z{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } -= $id_fangxiang{$id}{"1"};
				if( $count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) } == 0 ){
					delete $count{$chr}{$pos+$i}{$primer}{ substr($seq1,$i,1) }; }
				if(not exists$count{$chr}{$pos+$i}{$primer}{ substr($seq2,$i,1) }){
					print "$chr,$pos,$len,$i,$primer:",substr($seq2,$i,1),"\tnotadd\n"; next;}
				$count{$chr}{$pos+$i}{$primer}{ substr($seq2,$i,1) }--;
				$count_z{$chr}{$pos+$i}{$primer}{ substr($seq2,$i,1) } -= $id_fangxiang{$id}{"2"};
				if( $count{$chr}{$pos+$i}{$primer}{ substr($seq2,$i,1) }==0 ){
					delete $count{$chr}{$pos+$i}{$primer}{ substr($seq2,$i,1) };}
			}
		}
	}
}

sub print_file_header{
	my ($bam) = @_ ;
	print OUT "##fileformat=VCFv4.2\n",
		"##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
		"##samtoolsVersion=1.5+htslib-1.5\n",
		"##reference=file:///data/bioit/biodata/duyp/bin/hg19/hg19.fasta\n",
		"##contig=<ID=chr1,length=249250621>\n",
		"##contig=<ID=chr2,length=243199373>\n",
		"##contig=<ID=chr3,length=198022430>\n",
		"##contig=<ID=chr4,length=191154276>\n",
		"##contig=<ID=chr5,length=180915260>\n",
		"##contig=<ID=chr6,length=171115067>\n",
		"##contig=<ID=chr7,length=159138663>\n",
		"##contig=<ID=chr8,length=146364022>\n",
		"##contig=<ID=chr9,length=141213431>\n",
		"##contig=<ID=chr10,length=135534747>\n",
		"##contig=<ID=chr11,length=135006516>\n",
		"##contig=<ID=chr12,length=133851895>\n",
		"##contig=<ID=chr13,length=115169878>\n",
		"##contig=<ID=chr14,length=107349540>\n",
		"##contig=<ID=chr15,length=102531392>\n",
		"##contig=<ID=chr16,length=90354753>\n",
		"##contig=<ID=chr17,length=81195210>\n",
		"##contig=<ID=chr18,length=78077248>\n",
		"##contig=<ID=chr19,length=59128983>\n",
		"##contig=<ID=chr20,length=63025520>\n",
		"##contig=<ID=chr21,length=48129895>\n",
		"##contig=<ID=chr22,length=51304566>\n",
		"##contig=<ID=chrX,length=155270560>\n",
		"##contig=<ID=chrY,length=59373566>\n",
		"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n",
		"##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n",
		"##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">\n",
		"##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">\n",
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n",
		"##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">\n",
		"##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">\n",
		"##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">\n",
		"##INFO=<ID=BQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias (bigger is better)\">\n",
		"##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)\">\n",
		"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">\n",
		"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">\n",
		"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">\n",
		"##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">\n",
		"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n",
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n",
		"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">\n",
		"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n",
		"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">\n",
		"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">\n",
		"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">\n",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$bam\n";
}
my %ref_insert_genome;
my %ref_delete_genome;
sub check_indel{
	my($chr,$pos,$indel,$del_or_insert)=@_;	
	my($ref,$alt);
	if(exists $ref_insert_genome{$chr}{$pos-1}){
		my @ref_del = @{$ref_insert_genome{$chr}{$pos-1}};
		if( length($ref_del[1]) >= length($indel) ){ 
			$ref = $ref_del[0];
			$alt = join '',($ref_del[0],$indel); 
			return($ref,$alt);
		}else{}
	}elsif(exists $ref_delete_genome{$chr}{$pos-1-1}){
		my @ref_del = @{$ref_delete_genome{$chr}{$pos-1-1}};
		if( length($ref_del[1]) >= length($indel) ){ #print "@ref_del\t$indel\n";
			$ref = join '',@ref_del[0..1]; $alt = $ref;
			$alt =~ s/$indel//;
			return($ref,$alt);
		}else{print length($ref_del[1]),"\t<\t",length($indel),"\t$chr,$pos,$indel,$del_or_insert,$ref_del[1]\n";}
	}
	
	my $single_rep = &get_single_rep($indel);
	my $end = $pos + 40;
	my $sta = $pos - 1;
	if( $del_or_insert eq "delete" ){
		$sta = $sta - 1;
	}
	my $seq_ref_del = `samtools faidx $RefGenome $chr:$sta-$end`;
	$seq_ref_del  =~ tr/ATGCatgc/ATGCATGC/;
	my( $ref_id, $seq_ref_1, $seq_ref_2 ) = split( /\n/,$seq_ref_del);
	if($seq_ref_1 =~ /^(?<before>[ATGCN])(?<rep>($single_rep)+)(?<after>[ATGCN])/ and $del_or_insert eq "delete" ){ 
		my($before,$rep,$after) = ($+{before},$+{rep},$+{after});
		$ref = join '',($before,$rep); $alt = $ref;
		$alt =~ s/$indel//; 
		$ref_delete_genome{$chr}{$pos-1-1} = [$before,$rep,$after,$indel];
		# print "[$before,$rep,$after,$indel]\t$indel\t$single_rep\n";
	}elsif($seq_ref_1 =~ /^(?<before>[ATGCN])(?<rep>($single_rep)+)(?<after>[ATGCN])/ and $del_or_insert eq "insert" ){
		my($before,$rep,$after) = ($+{before},$+{rep},$+{after});
		$ref = join '',($before);
		$alt = join '',($before,$indel);
		$ref_insert_genome{$chr}{$pos-1} = [$before,"",$after,$indel];
	}else{
		my($before,$rep,$after) = (substr($seq_ref_1,0,1),$indel,substr($seq_ref_1,1,1));
		$ref = join '',($before);
		$alt = join '',($before,$indel);
		$ref_insert_genome{$chr}{$pos-1} = [$before,"",$after,$indel];
	}
	return($ref,$alt);
}

sub get_single_rep{
	my ($del,$single_rep) = @_;
	my @line = split //,$del;
	foreach my $i(0..($#line-1)){
		my $lll = join '',@line[0..$i];
		if( $del =~ /^($lll)+$/ ){
			$single_rep = $lll;
		}
	}
	if(defined $single_rep ){
		if( $single_rep ne $del){
			my $single_rep1 = &get_single_rep($single_rep);
		}else{return $single_rep;}
	}else{return $del;}
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0 (count base depth without overlap!)
Version: $version
Contact: Wu Guizhi<guizhi.wu\@genetalks.com> 

Usage:
  Options:
  -s  <sampleid>         sampleid, optional, default Sample
  -i  <file>             Input bam file, forced
  -o  <outfile>          Output file, forced
  -e  <endborderN>       Removed endborderN base, default 2 
  -r1 <read1only>        Used read1 only,[1:yes] 
  -d  <outdir>           Outdir, default ./
  -r  <regionfile>       Output region ,optional 
  -rp <regionprimerfile> Output region and output primer of each position ,optional
  -h                     Help

USAGE
	print $usage;
	exit;
}

