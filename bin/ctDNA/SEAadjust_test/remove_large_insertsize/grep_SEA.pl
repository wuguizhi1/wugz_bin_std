open(IN,"$ARGV[0]/$ARGV[1].fusion.final.txt");
while(<IN>){
	chomp;
	my$lines=$_;
	if($lines =~ /^chr16,215395,A/){
		my@line=split/\t/,$lines;
		my$ratio=$line[3]/$line[4];
		print "$ARGV[1]\tSEA\t$line[3]\t$line[4]\t$ratio\n";
	}
}
close IN;

