$dir_annovar_db="/data/bioit/biodata/database/humandb/";

sub read_data_sheet{
	my ($f, $adata)=@_;
	open (S, $f) or die $!;
	while (<S>) {
		chomp;
		next if(/^$/ || /^#/ || /^SampleName/);
		my ($sample, $index, $barInfo, $bed, $primerfile, $fregion_primer, $fusion_pos, $variant_target, $variant_before, @sample_info)=split /\t/, $_;
		$adata->{$sample}->{"barGroup"}=$barInfo;
		$adata->{$sample}->{"bedFile"}=$bed;
        $adata->{$sample}->{"primerFile"}=$primerfile;
        $adata->{$sample}->{"fregion_primer"}=$fregion_primer;
		my @variant_target = split /,/, $variant_target;
		if (defined $variant_target[0]) {
			$adata->{$sample}->{"variant_annovar"} = $variant_target[0];
		}
		if (defined $variant_target[1]) {
			$adata->{$sample}->{"variant_genotype"} = $variant_target[1];
		}
		$adata->{$sample}->{"index"}=$index;
		@{$adata->{$sample}->{"sample_info"}}=@sample_info;
		if ($variant_before ne "") {
			$adata->{$sample}->{"variant_before"}=$variant_before;
		}
		if ($fusion_pos ne "") {
			$adata->{$sample}->{"fusion_pos"}=$fusion_pos;
		}
	}
	close (S);
}

