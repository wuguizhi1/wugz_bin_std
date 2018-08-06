#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
####################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my($bamfile,$outdir,$outfile);
my($seq_error_rate,$min_simi,$minmismatch,$maxborder,$RefGenome,$print_type);
GetOptions(
		"help|?"	=>\&USAGE,
		"i:s"		=>\$bamfile,
		"o:s"   	=>\$outfile,
		"d:s"   	=>\$outdir,
		"ser:s"		=>\$seq_error_rate, # 根据测序错误率过滤一致性序列
		"mins:s"	=>\$min_simi,		# 相似度的底线值，低于则返回（0,0,0）
		"maxb:s"	=>\$maxborder,		# 参考序列边界截取长度
		"minmis:s"	=>\$minmismatch, 	# 根据错配数找边界
		"ref:s"		=>\$RefGenome,
		"p:s"		=>\$print_type,
) or &USAGE;
&USAGE unless defined($bamfile and $outfile);
if( not defined $seq_error_rate ){ $seq_error_rate = 0.005; }
if( not defined $min_simi ){	$min_simi = 0.75;	}
if( not defined $minmismatch ){	$minmismatch = 4;	}
if( not defined $maxborder ){	$maxborder = 80;	}
if( not defined $print_type ){  $print_type = 0; }
if( not defined $RefGenome ){   $RefGenome = "/data/bioit/biodata/duyp/bin/hg19/hg19.fasta"; }
if( not defined $outdir ){	$outdir = "./"; }
`mkdir $outdir` unless (-d $outdir);
$outdir=&AbsolutePath("dir",$outdir);

my $label_name = "WZ:Z:";
my( %genome_len, %DI_chr_pos, %finial_pos, %all_ref_count, %result_bam );
&SHOW_TIME("$bamfile\nstart"); ##
open(OUT,"|samtools view -bS - >$outfile");
&read_bam_head();
&read_I_D_of_bam();
&perfect_I_D_genome();
&read_I_D_of_bam_print();
close OUT;
&SHOW_TIME("end"); ##
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
	  	-i	<file>	Input bam file, forced
	  	-o	<file>	Output file, forced
	  	-d	<str>	Outdir, [./]
		-ser	<str>	seq_error_rate, default 0.005
		-mins	<str>	Min similarity, default 0.75
		-maxb	<int>	Max border, default 80
		-minmis	<int>	Min mismatch, default 4;
	  	-p	<int>	Print info in detail, default 0 [0:no] [1:yes]
	  	-h		Help

USAGE
	print $usage;
	exit;
}
sub AbsolutePath{		
	#获取指定目录或文件的决定路径
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
sub SHOW_TIME {
	#显示当时时间函数，参数内容是时间前得提示信息，为字符串
	my ($str)=@_;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	my $temp=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
	print "$str:\t[".$temp."]\n";
}
sub read_bam_head{
	open(I, "samtools view -H $bamfile|") or die $!;
	while(<I>){
		chomp;
		my $line = $_; 
		print OUT $_,"\n";
		if( $line =~ /SN:(.+)\tLN:(.+)$/ ){	#print "$1:$2\n";
			$genome_len{$1} = $2;
		}
	}
	close I;
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

### 按行读入文件，读取有DI的信息
sub read_I_D_of_bam{
	#&SHOW_TIME("Reading $bamfile"); ##
	open(I, "samtools view $bamfile|") or die $!;
	my $count_lines = 0;
	while(<I>){
		chomp; my$lines = $_;
		my ($id,$flag,$chr,$pos,$mapQ,$matchInfo,$li7,$li8,$li9,$seq,@line) = split(/\t/,$lines);
		my ($rnum, $is_reverse ,$is_unmap)=&explain_bam_flag($flag); ##		
		if ($is_unmap == 1 ){
			if( $print_type == 1 ){	print "$id,$flag,$chr,$pos\t This id is_unmap!\n";	}
		}else{
			$count_lines++;
			if( $matchInfo =~ m/[DI]/){
				&bianli_match_n_DI($matchInfo,$chr,$pos,$seq,"DI");
			}else{}
		}
	}
	close I;
	#&SHOW_TIME("Closing $count_lines"); ##
}
###	读DI 
sub bianli_match_n_DI{
	my($matchInfo,$chr,$pos,$seq,$M_pos)=@_;
	#print "$matchInfo,$chr,$pos,$seq,$M_pos\n";
	my($match_n_,$match_str_) = &sep_match_info($matchInfo);
	my @match_n = @$match_n_;
	my @match_str = @$match_str_;
	my ($seq_l,$seq_r)=(0,0);
	my($ppp,$next_str,$next_pos,$next_seq,$index);
	foreach my $index(0..$#match_n){
		if( $match_str[$index] eq "D" ){
			if( $M_pos eq "DI" ){
				$DI_chr_pos{$chr}{$seq_r + $pos - 1}{$match_n[$index]."D"}=1;
			}# 靠左的位置
			$seq_r += $match_n[$index];
			next;
		}
		$seq_l += $match_n[$index];
		if( $match_str[$index] eq "M" ){	$seq_r += $match_n[$index];	}
		if( $match_str[$index] eq "I" ){
			if( $index==0 or $index==$#match_n ){ next; }
			if( $M_pos eq "DI" ){ 
				$DI_chr_pos{$chr}{$seq_r + $pos - 1}{substr($seq,$seq_l-$match_n[$index],$match_n[$index])."I"}=1; 
			}
		}
	}
	return($pos,$pos+$seq_r-1);
}
### 每个DI单独存储，获取存储ref相关序列，并针对该DI所在位置，更改bam
sub perfect_I_D_genome{
	foreach my$chr(sort keys%DI_chr_pos){
		my$DI_chr = $DI_chr_pos{$chr};
		foreach my$pos(sort{$a<=>$b} keys%$DI_chr){
			my$DI_chr_ = $DI_chr_pos{$chr}{$pos};
			foreach my$ATGC_DI(sort keys %$DI_chr_){
				my $DI = substr($ATGC_DI,-1);
				my $alt = substr($ATGC_DI,0,length($ATGC_DI)-1); 
				if( $ATGC_DI =~ /D/ ){
					my ($start,$end) = ($pos+1,$pos+$alt);
					$alt = &get_samtools_faidx($chr,$start,$end);
				}
				my $count = $DI_chr_pos{$chr}{$pos}{$ATGC_DI};
				delete $DI_chr_pos{$chr}{$pos}{$ATGC_DI};
				$DI_chr_pos{$chr}{$pos}{$alt.$DI} = $count;
				&extend_refseq($chr,$pos,$DI,$alt);
			}
		}
	}
}
sub extend_refseq{
	my($chr,$pos,$DI,$alt)=@_;
	my( $start, $end ) = ($pos - $maxborder + 1, $pos + $maxborder);
	if( $pos - $maxborder +1 <= 0 ){	$start = 1; }
	if( $pos + $maxborder > $genome_len{$chr} ){		$end = $genome_len{$chr};	}

	my $seq_maxborder = &get_samtools_faidx($chr,$start,$end);
	chomp($seq_maxborder);
	my $single_rep = &get_single_rep($alt);		
	#print "$seq_maxborder\n$chr,$pos,$DI,$alt\t$single_rep\n";#
	my $before = substr($seq_maxborder, 0, $maxborder);
	my $repafter = substr($seq_maxborder, $maxborder);
	my ($rep, $after, $rep_);
	my $alt_len = 0;
	if($repafter =~ /^(?<rep>($single_rep)+)(?<after>[ATGCN]+)$/ ){
		($rep,$after) = ($+{rep},$+{after});
		if( $DI eq "D" ){ 
			$rep_ = $rep; $rep_ =~ s/$alt//; 
			$alt_len = length($alt);
		}
		if( $DI eq "I" ){ $rep_ = $alt.$rep;    }
	}else{
		($rep,$after) = ("",$repafter);
		$rep_ = $alt.$rep;
	}

	my $n_b = &extend_n_each($before.$rep_,$before.$rep,"B",$minmismatch);	#minmismatch bp差异
	my $n_a = &extend_n_each($rep.$after,$rep_.$after,"A",$minmismatch);
	my ($posl, $posr, $posrr) = ($pos, $pos+$alt_len+1, $pos+length($rep)+1);
	if( $print_type == 1 ){ print "$chr,$pos,$DI,$alt\n"; }
	my $each = "$chr,$start,$n_b,$posl,$posr,$posrr,$n_a,$end,$before,$rep_,$rep,$after,$DI,$alt";
	if( $print_type == 1 ){ print "$each\n";	}
	&cut_bam_for_ID($each);
}
### 指定位置，截取参考基因组序列
sub get_samtools_faidx{
	my($chr,$sta,$end)=@_;
	my $genome_ref = `samtools faidx $RefGenome $chr:$sta-$end`;
	$genome_ref =~ tr/ATGCatgc/ATGCATGC/;
	my( $ref_id, @seq_ref_1 ) = split( /\n/,$genome_ref );
	my $ref_seq12 = join "",@seq_ref_1;
	return($ref_seq12);
}
### 迭代得到指定序列的最小重复单元
sub get_single_rep{
	my ($alt,$single_rep) = @_;
	my @line = split //,$alt;
	foreach my $i(0..($#line-1)){
		my $smaller_alt = join "",@line[0..$i];
		if( $alt =~ /^($smaller_alt)+$/ ){
			$single_rep = $smaller_alt;
		}
	}
	if(defined $single_rep ){
		if( $single_rep ne $alt){
			my $single_rep1 = &get_single_rep($single_rep);
			return $single_rep1;
		}else{
			return $single_rep;
		}
	}else{
		return $alt;
	}
}
sub extend_n_each{
	my ($seq0,$seq1,$typ,$mis_n) = @_;
	if($typ eq "B"){
		$seq1 = reverse $seq1;
		$seq0 = reverse $seq0;
	}
	my @unit0 = split//,$seq0;
	my @unit1 = split//,$seq1;
	my $mis_n_new=0;
	foreach my$index(0..$#unit0){
		if( $unit0[$index] ne $unit1[$index] ){
			$mis_n_new++;
			if( $mis_n_new == $mis_n ){ 
				return($index+1); 
			}
		}
	}
}




sub cut_bam_for_ID{
	my( $each ) = @_;
	my($chr,$start,$n_b,$posl,$posr,$posrr,$n_a,$end,$before,$rep_,$rep,$after,$DI,$alt) = split/,/,$each; 
	my $len_D = 0;
	if( $DI eq "D" ){ $len_D = length($alt); }
	my $len_plus = length($rep);
	if( $DI eq "D" ){	$len_plus = $len_plus - length($alt); }
	my($left, $right) = ($posl-($n_b-$len_plus)+1, $posrr+($n_a-$len_plus)-1);
	my $genomeref = substr($before, -($n_b-$len_plus)).$rep.substr( $after,0,$n_a-$len_plus);
	my $genomealt = substr($before, -($n_b-$len_plus)).$rep_.substr($after,0,$n_a-$len_plus);
	$all_ref_count{$chr}{$posl}{$DI}{$alt}{"ref"}{$genomeref} = 0;
	$all_ref_count{$chr}{$posl}{$DI}{$alt}{"alt"}{$genomealt} = 0;
	my($mis_lenl,$max_lenl) = &similarity_mis($genomealt,$genomeref,"A");
	my($mis_lenr,$max_lenr) = &similarity_mis($genomealt,$genomeref,"B");
	$finial_pos{$posl.$DI.$alt}  = $each.",$mis_lenl,$max_lenl";
	$finial_pos{$posrr.$DI.$alt} = $each.",$mis_lenr,$max_lenr";
	my %store_bam;
	open(I, "samtools view $bamfile $chr:$posl\-$posrr|") or die $!;
	while(<I>){
		chomp; my$lines = $_;
		my ($id,$flag,$chr_,$pos,$mapQ,$matchInfo,$li7,$li8,$li9,$seq,@line) = split(/\t/,$lines);
		my($re,$se)=&before_bianli_match_n($matchInfo,$chr,$pos,$seq,$left,$right,$genomeref,$genomealt,$posl,$posrr,$DI,$alt,"A");
		if( $re ne "0" ){	$store_bam{$id}{$flag}{$re} = $se;	}
	}
	close I;
	my (@seq_ref,@seq_alt);
	my $ref_alt = $all_ref_count{$chr}{$posl}{$DI}{$alt};
	foreach my$alt_2(sort keys %$ref_alt){
		my $ref_count = $all_ref_count{$chr}{$posl}{$DI}{$alt}{$alt_2};
		my ($count_n,$count_first) = (0,0);
		foreach my$alt_each(sort{$ref_count->{$b}<=>$ref_count->{$a}} keys %$ref_count){ $count_n++;
			my $count = $all_ref_count{$chr}{$posl}{$DI}{$alt}{$alt_2}{$alt_each};
			if( $count_n == 1 ){ $count_first = $count; }
			if( $count < $seq_error_rate*$count_first and $count != 0 ){
				delete $all_ref_count{$chr}{$posl}{$DI}{$alt}{$alt_2}{$alt_each};
			}
			if( $print_type == 1 ){
print "COUNT:$chr,$posl,$DI,$alt\t$alt_2\t$alt_each\t$count\t$genomealt,$genomeref\t$mis_lenl,$max_lenl\t$mis_lenr,$max_lenr\n";
			}
		}
		if( $alt_2 eq "ref" ){ @seq_ref = sort{$ref_count->{$a}<=>$ref_count->{$b}}keys%$ref_count; }
		if( $alt_2 eq "alt" ){ @seq_alt = sort{$ref_count->{$a}<=>$ref_count->{$b}}keys%$ref_count; }
	}

	foreach my$id(sort keys%store_bam){
		my$store_bam_id = $store_bam{$id};
		foreach my$flag(sort keys%$store_bam_id){
			my ($seq_border,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave);
			if( exists $store_bam_id->{$flag}->{"A"} ){	
				$seq_border = $store_bam{$id}{$flag}{"A"};
				if( $print_type == 1 ){ print "$id,$flag,A,$mis_lenl,$max_lenl\t$seq_border\t#@seq_ref#@seq_alt#\t$DI,$alt\n";}	
				if( length($seq_border) <= $n_b-length($rep)+$len_D ){			next;
				}elsif( length($seq_border) <= $mis_lenl ){			
					($simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave)=(1,0,1,0,1,0,1,0,"");
				}else{
					($simi1,$mis1,$simi2,$mis2,$seq_leave)=&check_each_bam($seq_border,\@seq_ref,$mis_lenl,$max_lenl,"A");
					($simi3,$mis3,$simi4,$mis4)=&check_each_bam($seq_border,\@seq_alt,$mis_lenl,$max_lenl,"A");
				}
				if( $print_type == 1 ){ print "\t$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave\n";   }
				&check_each_simi($simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$id,$flag,"A",$seq_leave,
							$seq_border,$posl,$DI,$alt,$mis_lenl);
			}elsif( exists $store_bam_id->{$flag}->{"B"} ){	
				$seq_border = $store_bam{$id}{$flag}{"B"};	
				if( $print_type == 1 ){ print "$id,$flag,B,$mis_lenr,$max_lenr\t$seq_border\t#@seq_ref#@seq_alt#\t$DI,$alt\n"; }
				if( length($seq_border) <= $n_a-length($rep)+$len_D ){			next;
				}elsif( length($seq_border) <= $mis_lenr ){		
					($simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave)=(1,0,1,0,1,0,1,0,"");
				}else{
					($simi1,$mis1,$simi2,$mis2,$seq_leave)=&check_each_bam($seq_border,\@seq_ref,$mis_lenr,$max_lenr,"B");
					($simi3,$mis3,$simi4,$mis4)=&check_each_bam($seq_border,\@seq_alt,$mis_lenr,$max_lenr,"B");
				}
				if( $print_type == 1 ){ print "\t$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave\n";	}
				&check_each_simi($simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$id,$flag,"B",$seq_leave,
							$seq_border,$posrr,$DI,$alt,$mis_lenr);
			}
		}
	}
}
sub check_each_simi{
	my($simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$id,$flag,$goon,$seq_leave,$seq_border,$pos,$DI,$alt)=@_;
	my $re = "";
	if( $mis2 - $mis4 >= $minmismatch-1 or $simi4 - $simi2 == 1 ){
		$re = "alt";
	}elsif( $mis4 - $mis2 >= $minmismatch-1 or $simi2 - $simi4 == 1 ){
		$re = "ref";
	}else{
		$re = "no";
	}
	my @re_new = ($re,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave,$seq_border,$pos,$DI,$alt);
	if( $print_type == 1 ){	print "\t\t@re_new\n";	}
	if( exists $result_bam{$id}{$flag}{$goon} ){
		my@result = @{$result_bam{$id}{$flag}{$goon}};
		if( $result[0] eq "no" and $re ne "no" ){	# old  no
			@{$result_bam{$id}{$flag}{$goon}} =
				($re,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave,$seq_border,$pos,$DI,$alt);
		}elsif( $result[0] ne "no" and $re eq "no" ){ # new no
		}elsif( $result[0] eq "no" and $re eq "no" ){	# all no
			if( $goon eq "A" and $result[11] < $pos ){ 			# 取前一判断
			}elsif( $goon eq "B" and $result[11] > $pos  ){  	# 取前一判断
			}elsif( $result[11] == $pos ){						# 不用比较
			}else{  
				if( $print_type == 1 ){	print "get this\t";	}
				@{$result_bam{$id}{$flag}{$goon}} = 
					($re,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave,$seq_border,$pos,$DI,$alt);
			}
		}elsif( $result[0] ne "no" and $re ne "no" ){	# all not no
			if( $goon eq "A" and $result[11] < $pos ){			# 取前一判断
			}elsif( $goon eq "B" and $result[11] > $pos  ){		# 取前一判断
			}elsif( $result[11] == $pos ){						# 取高分,old高分
				if( $result[0] eq "alt" and $re eq "alt" and $result[5]+$result[7]>=$simi3+$simi4){	
				}elsif( $result[0] eq "alt" and $re eq "ref" and $result[5]+$result[7]>=$simi1+$simi2){
				}elsif( $result[0] eq "ref" and $re eq "ref" and $result[1]+$result[3]>=$simi1+$simi2){
				}else{
					if( $print_type == 1 ){ print "get this\t"; }
					@{$result_bam{$id}{$flag}{$goon}} =
						($re,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave,$seq_border,$pos,$DI,$alt);
				}
			}else{
				if( $print_type == 1 ){ print "get this\t"; }
				@{$result_bam{$id}{$flag}{$goon}} =
					($re,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave,$seq_border,$pos,$DI,$alt);
			}
			if( $print_type == 1 ){ print "\t what:$re,$pos,$DI,$alt\t$result[0],$result[11],$result[12],$result[13]\n";	}
		}
		if( $print_type == 1 ){
			if( $result[0] ne $re ){
				print "before notsame:\t@result\n";
			}else{
				print "before same:\t@result\n";
			}
		}
	}else{
		@{$result_bam{$id}{$flag}{$goon}} = 
				($re,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave,$seq_border,$pos,$DI,$alt);
	}
}
sub check_each_bam{
	my($seq_border,$seq_alt,$mis_len,$max_len,$goon) = @_;
	my($simi1,$mis1,$simi2,$mis2);
	my $seq_leave = "";
	foreach my$seq_each(@$seq_alt){
		my ($seq_1, $seq_1_, $seq_2, $seq_2_);
		if( $goon eq "A" ){
			$seq_1 = substr($seq_border, 0, $mis_len);
			$seq_1_= substr($seq_each, 0, $mis_len);
			$seq_2 = substr($seq_border, $mis_len, $max_len-$mis_len);
			$seq_2_= substr($seq_each, $mis_len, $max_len-$mis_len);
		}else{
			$seq_1 = substr($seq_border, -$mis_len);
			$seq_1_= substr($seq_each, -$mis_len);
			$seq_2 = substr($seq_border, -$max_len);
			$seq_2_= substr($seq_each, -$max_len);
			if( length($seq_2) < $max_len ){
				$seq_2 = substr($seq_2, 0, length($seq_2)-$mis_len);
			}else{
				$seq_2 = substr($seq_2, 0, $max_len-$mis_len);
			}
			$seq_2_= substr($seq_2_,0, $max_len-$mis_len);
		}
		my($simi1_, $mis1_)=&similarity_simi($seq_1, $seq_1_, $goon);
		my($simi2_, $mis2_)=&similarity_simi($seq_2, $seq_2_, $goon);
		if( not defined $mis1 ){ ($simi1,$mis1,$simi2,$mis2) = ($simi1_, $mis1_, $simi2_, $mis2_); }
		elsif( $mis1_ < $mis1 ){ ($simi1,$mis1) = ($simi1_, $mis1_); }
		elsif( $mis2_ < $mis2 ){ ($simi2,$mis2) = ($simi2_, $mis2_); }
		if( $print_type == 1 ){	print "\t\t$simi1,$mis1,$simi2,$mis2\t$seq_1, $seq_1_, $goon, $seq_2, $seq_2_\n";	}
	}
	if( length($seq_border) > $max_len and $goon eq "A"  ){ 
		$seq_leave = substr($seq_border,-(length($seq_border)-$max_len));
	}elsif( length($seq_border) > $max_len and $goon eq "B"  ){
		$seq_leave = substr($seq_border, 0, (length($seq_border)-$max_len));
	}
	return($simi1,$mis1,$simi2,$mis2,$seq_leave);
}
sub before_bianli_match_n{
	my ($matchInfo,$chr,$pos,$seq,$left,$right,$genomeref,$genomealt,$posl,$posr,$DI,$alt,$goon) = @_;
	my ($plenl,$prefl,$pnl,$pindexl) = &bianli_match_n($matchInfo,$chr,$pos,$seq,$posl,$alt);
	my ($plenr,$prefr,$pnr,$pindexr) = &bianli_match_n($matchInfo,$chr,$pos,$seq,$posr,$alt);
	my ($lenl,$refl,$nl,$indexl,$nextlstr,$nextln,$nextlseq) = &bianli_match_n($matchInfo,$chr,$pos,$seq,$left,$alt);
	my ($lenr,$refr,$nr,$indexr,$nextrstr,$nextrn,$nextrseq) = &bianli_match_n($matchInfo,$chr,$pos,$seq,$right,$alt);
	if( defined $indexl and defined $indexr ){
		if( $indexr - $indexl == 0 ){
			my $ref_start = $left - ($pos + $refl);
			my $ref_end = $right - ($pos + $refl);
			my $ref_each = substr($seq, $ref_start+$lenl, $ref_end-$ref_start+1);
			if( $genomeref ne $ref_each ){
				my $simi_len = &similarity_mis($genomeref,$ref_each,$goon);
				if( $print_type == 1 ){
					if( length($genomeref) eq $simi_len ){ print "True\n"; }
					print "ref:$left,$right\n# $genomeref\n# $ref_each\t$lenl+$ref_start,$lenr+$ref_end\n";
					print "$lenl,$refl,$nl,$indexl,$nextlstr,$nextlseq\n$lenr,$refr,$nr,$indexr,$nextrstr,$nextrseq\n"; 
				}
			}	
			$all_ref_count{$chr}{$posl}{$DI}{$alt}{"ref"}{$ref_each}++;
		}elsif( $indexr - $indexl == 2  ){
			my $ref_start = $left - ($pos + $refl);
			my $ref_end = $right - ($pos + $refr);
			my $alt_each = substr($seq, $ref_start+$lenl, ($ref_end+$lenr) - ($ref_start+$lenl)+1);
			if( $genomealt ne $alt_each ){ 
				my $simi_len = &similarity_mis($genomealt,$alt_each,$goon);
				if( $print_type == 1 ){
					if( length($genomealt) eq $simi_len ){ print "True\t"; }
					print "alt:$left,$right\n## $genomealt\n## $alt_each\t$lenl+$ref_start,$lenr+$ref_end\n"; 
					print "$lenl,$refl,$nl,$indexl,$nextlstr,$nextlseq\n$lenr,$refr,$nr,$indexr,$nextrstr,$nextrseq\n";
					if( $nextlseq ne $alt or $nextlstr ne $DI ){
						print "####$nextlstr eq $DI \t $nextlseq eq $alt\n";
					}
				}
			}	
			if( $DI eq "D" ){
				if( $nextlseq eq $alt and $nextlstr eq $DI ){
					$all_ref_count{$chr}{$posl}{$DI}{$alt}{"alt"}{$alt_each}++;
				}
			}else{
				if( $nextlseq eq $alt and $nextlstr eq $DI ){
					$all_ref_count{$chr}{$posl}{$DI}{$alt}{"alt"}{$alt_each}++;
				}
			}
		}else{
			#die "error:not ref not alt!\t$matchInfo,$chr,$pos\t$left,$right,$posl,$posr\n";
		}
		return( 0 );
	}elsif( not defined $pindexr and defined $indexl ){
		if( not defined $refl ){ die "error:$matchInfo,$chr,$pos\t$left,$right,$posl,$posr\n"; }
		my $ref_start = $left - ($pos + $refl) + $lenl;
		my $seq_border = substr( $seq,-(length($seq)-$ref_start) );
		return( "A",$seq_border );
	}elsif( not defined $pindexl and defined $indexr ){
		my $ref_end = $right - ($pos + $refr) + $lenr + 1 ;
		my $seq_border = substr( $seq, 0, $ref_end );
		return( "B",$seq_border );
	}
}
sub similarity_simi{
	my($seq0, $seq1, $goon)=@_;
	my $lenm = length($seq0);
	if( length($seq0) > length($seq1) ){ $lenm = length($seq1); }
	if( $lenm == 0 ){
		return(0,0);
	}
	if($goon eq "B"){
		$seq1 = reverse $seq1;
		$seq0 = reverse $seq0;
	}
	$seq1 = substr($seq1,0,$lenm);
	$seq0 = substr($seq0,0,$lenm);
	my @unit0 = split //, $seq0;
	my @unit1 = split //, $seq1;

	my ($mismatch, $simi) = (0,0);
	foreach my$index(0..$#unit0){
		if( $unit0[$index] ne $unit1[$index] and $unit0[$index] ne "N" ){ # seq_border 非 N 才记为错配！
			$mismatch++;
		}
	}
	return(($lenm-$mismatch)/$lenm,$mismatch);
}
sub similarity_mis{
	my($seq0,$seq1,$goon)=@_;
	my $lenm = length($seq0);
	if( length($seq0) > length($seq1) ){	$lenm = length($seq1);	}
	if($goon eq "B"){
		$seq1 = reverse $seq1;
		$seq0 = reverse $seq0;
	}
	$seq1 = substr($seq1,0,$lenm);
	$seq0 = substr($seq0,0,$lenm);
	my @unit0 = split //, $seq0;
	my @unit1 = split //, $seq1;
	foreach my$index(0..$#unit0){
		if( $unit0[$index] ne $unit1[$index] ){
			return($index,$lenm);
		}
	}
	return($lenm,$lenm);
}
### 读指定位点在bam中的信息
### 读bam左右M端点
sub bianli_match_n{
	my($matchInfo,$chr,$pos,$seq,$M_pos,$alt)=@_;
	my($match_n_,$match_str_) = &sep_match_info($matchInfo);
	my @match_n = @$match_n_;
	my @match_str = @$match_str_;
	my ($seq_l,$seq_r,$seq_D,$seq_I,$seq_S,$seq_M)=(0,0,0,0,0,0);
	my($next_str,$next_seq,$next_n);
	foreach my $index(0..$#match_n){
		if( $match_str[$index] eq "D" ){
			if( $M_pos>=$seq_r+$pos and $M_pos<$match_n[$index]+$pos+$seq_r ){
				my($next_str,$next_seq,$next_n)=&cigar_index_next(\@match_n,\@match_str,$index,$alt,$seq,$seq_l,$seq_r+$pos,$chr);
				return( $seq_l,$seq_r,$match_n[$index],$index,$next_str,$next_n,$next_seq);
			}
			$seq_r += $match_n[$index];
			$seq_D += $match_n[$index];
			next;
		}
		if( $match_str[$index] eq "M" ){
			if( $M_pos>=$seq_r+$pos and $M_pos<$match_n[$index]+$pos+$seq_r ){
				my($next_str,$next_seq,$next_n)=&cigar_index_next(\@match_n,\@match_str,$index,$alt,$seq,$seq_l,$seq_r+$pos,$chr);
				return( $seq_l,$seq_r,$match_n[$index],$index,$next_str,$next_n,$next_seq);
			}
			$seq_r += $match_n[$index]; 
			$seq_l += $match_n[$index];
			$seq_M += $match_n[$index];
		}elsif( $match_str[$index] eq "I" ){
			$seq_l += $match_n[$index];
			$seq_I += $match_n[$index];
		}else{
			$seq_l += $match_n[$index];
			$seq_S += $match_n[$index];
		}
	}
}
sub cigar_index_next{
	my($match_n_,$match_str_,$index,$alt,$seq,$seq_l,$seq_r,$chr) = @_;
	my @match_n = @$match_n_;
	my @match_str = @$match_str_;
	my($next_str,$next_seq,$next_n);
	if( $#match_n > $index ){
		$next_str = $match_str[$index+1];
		$next_n = $match_n[$index+1];
		$next_seq = $match_n[$index+1].$match_str[$index+1];
		if( $next_str eq "I" ){ 
			$next_seq = substr($seq,$seq_l+$match_n[$index],$match_n[$index+1]); 
		}elsif( $next_str eq "D"){
			$next_seq = &get_samtools_faidx($chr,$seq_r+$match_n[$index],$seq_r+$match_n[$index]+$next_n-1);
		}
	}else{
		($next_str,$next_seq,$next_n) = ("","",0);
	}
	return($next_str,$next_seq,$next_n);
}
sub sep_match_info{
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
sub read_I_D_of_bam_print{
	open(I, "samtools view $bamfile|") or die $!;
	while(<I>){
		chomp; my$lines = $_;
		my ($id,$flag,@line) = split(/\t/,$lines);
		if( exists $result_bam{$id}{$flag} ){
			my $goon_each = $result_bam{$id}{$flag};
			foreach my$goon(sort keys %$goon_each){
				if( $print_type == 1 ){	print "\ngoon:$goon\t@{$goon_each->{$goon}}\n$lines\n";	}
				$lines = &change_I_D($id,$flag,$goon,$lines);
			}
			print OUT $lines,"\n";
		}else{
			print OUT $lines,"\n";
		}
	}
	close I;
}
sub change_I_D{
	my ($id,$flag,$goon,$lines)=@_;
	my ($lines_new, $matchInfo) = &change_cigar($id,$flag,$goon,$lines);
	return($lines_new);
}
sub change_cigar{
	my ($id,$flag,$goon,$lines)=@_;
	my $add_label_or_not;
	my @line = split(/\t/,$lines);
	my ($pos , $matchInfo) = ($line[3], $line[5]);
	my ($result,$simi1,$mis1,$simi2,$mis2,$simi3,$mis3,$simi4,$mis4,$seq_leave,$seq_border,$posM,$DI,$alt)=
			@{$result_bam{$id}{$flag}{$goon}};
	my $each = $finial_pos{$posM.$DI.$alt};
	my ($chr,$start,$n_b,$posl,$posr,$posrr,$n_a,$end,$before,$rep_,$rep,$after,$DI_,$alt_,$mis_len,$max_len) = split/,/,$each;
	my ($cigar1,$cigar2) = &sep_cigar($matchInfo,$posM,$pos); # 分割时posM 计入cigar1，拼接时需注意
	if( $cigar1 eq "0" ){ print $lines,"\n";return($lines,$matchInfo);  } # 分割时在D区域直接返回不处理
	my $len_MD_old1 = &len_cigar($cigar1,"MD");
	if( $print_type == 1 ){	print "$cigar1,$cigar2\n";	}
	my ($len_D, $len_I) = (0,0);
	if( $DI eq "D" ){ $len_D = length($alt); }
	if( $DI eq "I" ){ $len_I = length($alt); }
	my ($seq_border_A, $seq_border_B) = ("", "");
	if($n_b-length($rep)+$len_D <= length($seq_border) ){
		$seq_border_A = substr($seq_border, $n_b-length($rep)+$len_D );
	}
	if( (length($seq_border) - $n_a+length($rep)-$len_D) > 0 ){
		$seq_border_B = substr($seq_border, 0, (length($seq_border) - $n_a+length($rep)-$len_D) );
	}
	my ($A_leave_ref,$A_leave_alt,$B_leave_ref,$B_leave_alt) = ($rep.$after,$rep_.$after,$before.$rep,$before.$rep_);
	if( $print_type == 1 ){
		print "$before,$rep_,$rep,$after\n$A_leave_ref\n$A_leave_alt\n$seq_border_A\n$B_leave_ref\n$B_leave_alt\n$seq_border_B\n";
	}
	my $len_S1 = &len_cigar($cigar1,"MSI");
	my $len_S2 = &len_cigar($cigar2,"MSI");
	my $len_border = length($seq_border);
	my ($simi, $match, $matchs) = ("","","");
	# 需注意 $matchs 为零的情况！ 
	if( $result eq "no" and $goon eq "A" ){			$cigar2 = $len_S2."S";
	}elsif( $result eq "ref" and $goon eq "A" ){	
		($simi, $match, $matchs) = &get_similarity($goon, $seq_border_A, $A_leave_ref ,1);
		if( $len_S2 == $matchs  ){	
			$cigar2 = $len_S2."M";
		}elsif( $len_S2 > $matchs and $matchs > 0 ){ # $matchs>=2
			$cigar2 = ($matchs)."M".($len_S2-$matchs)."S";
		}elsif( $matchs == 0 ){
			$cigar2 = $len_S2."S";
		}else{ 
print "\t what:$line[0],$line[1]\trefA\tsimi$simi,$match,$matchs\t$len_S2,I$len_I, $seq_border_A, $A_leave_alt\n"; }
	}elsif( $result eq "alt" and $goon eq "A" ){    
		($simi, $match, $matchs) = &get_similarity($goon, $seq_border_A, $A_leave_alt ,1);
		if( $len_S2 <= $len_I ){	$cigar2 = $len_S2.$DI;
			if( $len_S2 < $len_I ){	$add_label_or_not = $label_name.length($alt).$DI."_".$alt; 	}
		}elsif( $len_S2 == $matchs ){
			$cigar2 = length($alt).$DI.($len_S2-$len_I)."M";
		}elsif( $len_S2 > $matchs and $matchs > 0 ){
			$cigar2 = length($alt).$DI.($matchs-$len_I)."M".($len_S2-$matchs)."S";
		}elsif( $matchs == 0 ){
			$cigar2 = $len_S2."S";
		}else{ 
print "\t what:$line[0],$line[1]\taltA\tsimi$simi,$match,$matchs\t$len_S2,I$len_I, $seq_border_A, $A_leave_alt\n"; }
	}elsif( $result eq "no" and $goon eq "B" ){		$cigar1 = ($len_S1-1)."S1M";
	}elsif( $result eq "ref" and $goon eq "B" ){	
		($simi, $match, $matchs) = &get_similarity($goon, $seq_border_B, $B_leave_ref ,1);
		if( $len_S1 == $matchs+1 ){	
			$cigar1 = $len_S1."M";
		}elsif( $len_S1 - $matchs - 1 > 0 ){	
			$cigar1 = ($len_S1-$matchs-1)."S".($matchs+1)."M";
		}else{ 
print "\t what:$line[0],$line[1]\trefB\tsimi$simi,$match,$matchs\t$len_S1,I$len_I, $seq_border_B, $B_leave_alt\n"; }
	}elsif( $result eq "alt" and $goon eq "B" ){
		($simi, $match, $matchs) = &get_similarity($goon, $seq_border_B, $B_leave_alt ,1);    
		if( $len_S1 <= $len_I + $posrr-$posr+1 ){	$cigar1 = ($len_S1-($posrr-$posr+1)).$DI.($posrr-$posr+1)."M";
			if( $len_S1 < $len_I + $posrr-$posr+1 ){ $add_label_or_not = $label_name.length($alt).$DI."_".$alt;   }
		}elsif( $len_S1 == $matchs + 1){
			$cigar1 = ($len_S1-$len_I-($posrr-$posr+1))."M".length($alt).$DI.($posrr-$posr+1)."M";
		}elsif( $len_S1 > $matchs + 1 ){
			$cigar1 = ($len_S1-($matchs)-1)."S".($matchs-$len_I-($posrr-$posr))."M".length($alt).$DI.($posrr-$posr+1)."M";
		}else{ 
print "\t what:$line[0],$line[1]\taltB\tsimi$simi,$match,$matchs\t$len_S1,I$len_I, $seq_border_B, $B_leave_alt\n"; }
	}
	if( $cigar2 =~ /-/ ){ $cigar2 = $len_S2."S";
	}elsif( $cigar1 =~ /-/ ){	$cigar1 = ($len_S1-1)."S1M";	}

	if( $print_type == 1 ){	print "$cigar1,$cigar2\tsimi:$simi, $match, $matchs\t$len_S2,$n_a,a,b,$n_b,$len_S1\n";	}
	my $cigar_new = &paste_cigar($cigar1,$cigar2);	
	my $len_MD_old = &len_cigar($matchInfo,"MD");
	my $len_MD_new = &len_cigar($cigar_new,"MD");
	my $len_MD_new1 = &len_cigar($cigar1,"MD");
	$line[3] = $line[3] + $len_MD_old1 - $len_MD_new1;
	if( $matchInfo ne $cigar_new and $print_type >= 1){	
		print "$chr\t$posM\t$goon\t$DI\t$alt\tcigar:$matchInfo\tnew:$cigar_new\t$line[0],$line[1]\n";		
	}elsif( $print_type == 1 ){	
		print "cigar:$matchInfo\told:$cigar_new\n"; 	
	}

	my $len_MSI_old = &len_cigar($matchInfo,"MSI");
	my $len_MSI_new = &len_cigar($cigar_new,"MSI");
	if( $len_MSI_old != $len_MSI_new ){	die "error:$line[0]\t$line[1]\t$matchInfo to $cigar_new?\n";	}
	if( $cigar_new =~ /-/ ){	
		print "\t what:$line[0],$line[1]\tcigar:$matchInfo\tnew:$cigar_new\thas - illegal cigar\n";	
	}elsif($cigar_new =~ /^0[MID]/ or $cigar_new =~ /[MID]0[MID]$/){
		print "\t what:$line[0],$line[1]\tcigar:$matchInfo\tnew:$cigar_new\thas 0 illegal cigar\n";
	}else{
		if( $matchInfo ne $cigar_new ){
			$line[5] = $cigar_new;
			$lines = join "\t",@line;
			$lines = &change_MD_to_no($lines);
		}
		if( defined $add_label_or_not ){	$lines = $lines."\t".$add_label_or_not; }
	}
	return($lines,$matchInfo);
}
sub change_MD_to_no{
	my ( $lines , $bef,$MD,$aft) = @_;
	if ( $lines =~ m/^(?<bef>.+)\tMD:Z:(?<MD>.+?)\t(?<aft>.+?)$/){
		($bef,$MD,$aft) = ($+{bef},$+{MD},$+{aft});
	}
	if( defined $MD ){
		$lines = $bef."\t".$aft;
	}
	return($lines);
}
### 2, 得到最长完全匹配数
### 1，不要求完全匹配，返回匹配度与匹配数；
### 0，连3不匹配时，提前返回匹配度与匹配数；
sub get_similarity_{
	my ($re , $seq, $seq0 ) = @_;
	if( length($seq)==0 ){ return(0,0); }
	my ($lenm) = &min_array(length($seq0), length($seq));
	my @unit0 = split //, $seq0;
	my @unit = split //, $seq;
	my ($match, $mismatch)=(0,0);
	for(my $i=0; $i<$lenm; $i++){
		if($unit[$i] eq $unit0[$i]){
			$match++;
		}else{
			$mismatch++;
			if( $re==2 ){ return(1,$match); }
		}
		if( $re==0 and $i>0 and $i<scalar(@unit)-2 ){ # 连续3个snp
			if($unit[$i] ne $unit0[$i] and $unit[$i+1] ne $unit0[$i+1] and $unit[$i+1] ne $unit0[$i+1] ){
				last;
			}
		}
	}

	my $simi = $match/($match+$mismatch);
	return ($simi, $match);
}
### 内部择优得到最佳匹配长度，并返回其对应的相似度，匹配数，截断序列长度
sub get_similarity{
	my ($ori, $seq, $seq0, $typ)=@_;
	my $lenm = length($seq0);
	if( length($seq0) > length($seq) ){ $lenm = length($seq); }
	if($ori eq "B"){
		$seq = reverse $seq;
		$seq0 = reverse $seq0;
	}
	$seq = substr($seq,0,$lenm);
	$seq0 = substr($seq0,0,$lenm);
	my%range;
	foreach my$i(0..length($seq)){
		my($simi1, $match1)=&get_similarity_(1,substr($seq,0,$i),substr($seq0,0,$i));
		my($simi2, $match2)=&get_similarity_(1,substr($seq,$i-length($seq)),substr($seq0,$i-length($seq)));
		if( $i == length($seq) ){ ($simi2, $match2) = (0,0); }
		@{$range{3*$simi1-$simi2}{$i}} = ($simi1, $match1, $simi2, $match2, substr($seq,0,$i),substr($seq0,0,$i),
				substr($seq,$i-length($seq)),substr($seq0,$i-length($seq)));
	}
=pob
	foreach my$cha(sort{$b<=>$a}  keys%range){
		my $range_max = $range{$cha};
		foreach my$i(sort{$a<=>$b}  keys%$range_max){
			print "\t$cha\t$i\t@{$range{$cha}{$i}}\n";
		}
	}
=cut
	my @range_p = keys%range;
	my($max, $max_n, $equal, $same1 ,$same2) = &max_array(@range_p);
	my $range_max = $range{$max};
	foreach my$i(sort{$a<=>$b}  keys%$range_max){
		my( $simi1, $match1, $simi2, $match2, $seq_, $seq0_ )=@{$range{$max}{$i}};
		if( $simi1 > $min_simi ){
			return ($simi1, $match1, $i);
		}else{
			return (0,0,0);
		}
	}
	return (0, 0, 0);
}
### 返回数组最大值，最大值位置，最大值是(0)否(1)唯一，第1个最大值位置，第2个最大值位置
sub max_array{
	my ($max,$max_n,$equal,$same1,$same2,$count);
	foreach (@_){
		my$a=$_;	$count++;
		if( $count==1 ){
			$max = $a;
			$max_n = $count;
			$equal = 0;
			$same1 = $count;
			$same2 = $count;
		}elsif( $a>$max ){
			$max = $a;
			$max_n = $count;
			$equal = 0;
			$same1 = $count;
			$same2 = $count;
		}elsif( $a==$max ){
			$equal = 1;
			$same2 = $count;
		}
	}
	return($max,$max_n,$equal,$same1,$same2);
}
### 返回数组最小值，最小值位置，最小值是(0)否(1)唯一，第1个最小值位置，第2个最小值位置
sub min_array{
	my ($min,$min_n,$equal,$same1,$same2,$count);
	foreach (@_){
		my$a=$_;	$count++;
		if( $count==1 ){
			$min = $a;
			$min_n = $count;
			$equal = 0;
			$same1 = $count;
			$same2 = $count;
		}elsif( $a<$min ){
			$min = $a;
			$min_n = $count;
			$equal = 0;
			$same1 = $count;
			$same2 = $count;
		}elsif( $a==$min ){
			$equal = 1;
			$same2 = $count;
		}
	}
	return($min,$min_n,$equal,$same1,$same2);
}
### 计数各个类型的长度
sub len_cigar{
	my($cigar,$typ)=@_;
	my ($match_n_,$match_str_) = &sep_match_info($cigar);
	my @match_n = @$match_n_;
	my @match_str = @$match_str_;
	my ($re,$re1) = (0,0);
	foreach my$index(0..$#match_n){
		if( $match_str[$index] =~ /[MSI]/ ){
			$re += $match_n[$index];	# 对应测序序列长度
		}
		if( $match_str[$index] =~ /[MD]/ ){
			$re1 += $match_n[$index];	# 对应参考序列长度
		}
	}
	if( $typ eq "MSI" ){
		return($re);
	}elsif( $typ eq "MD" ){
		return($re1);
	}
}
### 2个cigar的粘贴
sub paste_cigar{
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
sub sep_cigar{
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
				#print "in delete $matchInfo !\n";
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

