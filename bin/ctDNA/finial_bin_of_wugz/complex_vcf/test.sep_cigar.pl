my($matchInfo,$M_pos,$pos)=("66M",9,9);
my($c1, $c2) = &sep_cigar($matchInfo,$M_pos,$pos);
print "$matchInfo,$M_pos,$pos\t#$c1, #$c2\n";

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
