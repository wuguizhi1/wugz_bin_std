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
my ($fIn,$fkey,$outdir, $primer);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"p:s"=>\$primer,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my %end;
my %seq;
my %cigar;
my %end_start;
open(I, $fIn) or die $!;
while(<I>){
	chomp;
	next if(/^$/);
	my ($id, $flag, $chr, $pos, $qua, $cigar, undef, undef, undef, $seq)=split /\t/, $_;
	my ($pid)=$_=~/XP:Z:(\S+)/;
	next if($chr ne "chr16" || (defined $primer and $pid ne $primer));
	my ($is_read1, $strand) = &explain_strand_flag($flag);
	next if($is_read1==0 && $strand eq "+"); #
	my ($end_pos, $start_pos) = &get_end_pos($pos, $cigar, $strand);
	push @{$end{$id}}, $end_pos;
	$end_start{$id}{$end_pos} = $start_pos;
	if($is_read1){
		$seq{$id}=$seq;
		$cigar{$id}=$cigar;
	}
	print "$id\t$flag\t$pos\t$cigar\t$strand\t$end_pos\t$_\n";
}
close(I);

open(O,">$outdir/$fkey.size") or die $!;
foreach my $id (sort {$a cmp $b} keys %end){
	my @pos = @{$end{$id}};
	if(@pos ==1){
		print "one pos:",$id,"\n";
		next;
	}
	
	my @pos_sort = sort {$a <=> $b} @pos;
	my $len;
	my $type;
	if($seq{$id}=~/AGGTTCAC/){
		$type = "SEA";
	}elsif($seq{$id}=~/AGGTTCTA/){
		$type = "Ref";
	}else{
		print "unkown type:",$id,"\t",$seq{$id},"\n";
		next;
	}
	if($pos_sort[-1] >= 234700){
		$len = $pos_sort[-1] - 234700;
	}else{
		$len = $pos_sort[-1] - 215395;
	}
	
	
	if($len>1000){
		print "too long: $id\t$len\n";
		next;
	}
	if($len == 0){
		my ($s)=$cigar{$id}=~/(\d+)S/;
		$len = defined $s? $s: 0;
	}
	my $start_end = $end_start{$id};
	my @k = keys %$start_end;
	my @v = values %$start_end;
	print join("\t",$id, $type,$len),"\t#@k\t#@v\n";
	print O join("\t",$id, $type,$len),"\n";
}
close(O);



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub explain_strand_flag{
	my ($flag)=@_;
	my $flag_bin=sprintf("%b", $flag);
	my @flag_bin = split //, $flag_bin;
	my $is_reverse = $flag_bin[-5];
	my $strand = $is_reverse? "-": "+";
	my $is_read1 = $flag_bin[-7];
	return ($is_read1, $strand);
}
sub cigar_split{
	my($cigar)=@_;
	my @ucigar = split //, $cigar;
	my (@match_n,@match_str);
	my $nstr="";
	foreach my $i(0..$#ucigar){
		if($ucigar[$i] eq "M" || $ucigar[$i] eq "I" || $ucigar[$i] eq "D" || $ucigar[$i] eq "S"){
			push @match_str, $ucigar[$i];
			push @match_n, $nstr;
			$nstr = "";
		}else{
			$nstr .= $ucigar[$i];
		}
	}
	return(\@match_n,\@match_str);
}

sub get_end_pos{
	my ($pos,  $cigar, $strand)=@_;

	my ($cigar_n, $cigar_str)=&cigar_split($cigar); ## return the add of array
	my @cn = @{$cigar_n};
	my @cstr = @{$cigar_str};
	my $posg = $pos;
	my $posr = 0;
	my $sindex;
	my $bpos;

	for(my $i=0; $i<@cstr; $i++){
		if($cstr[$i] eq "S"){
			if($i==0){
				$posr+=$cn[$i];
			}else{
				$posg--;
			}
			$sindex = $i;
			$bpos = $posg;
		}elsif($cstr[$i] eq "M" ){
			$posg+=$cn[$i];
			$posr+=$cn[$i];
		}elsif($cstr[$i] eq "D" ){
			$posg+=$cn[$i];
		}elsif($cstr[$i] eq "I"){
			$posr+=$cn[$i];
		}
	}
	if( $strand eq "+"){
		return ($pos,$posg);
	}else{
		return ($posg,$pos);
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
  -i  <file>   Input file, forced
  -k  <str>	Key of output file, forced
  -p  <str>    primer id, optinal
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

