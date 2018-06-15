#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]);
while (<IN>) {
	chomp;
	@t = split(/\t/);
	$gene = "na";
	if($t[7] ne "."){
		@tt = split(/\|/, $t[7]);
		$gene = $tt[0];
	}
	unless (exists $hash{$t[3]}) {
		$hash{$t[3]} = [];
	}
	push (@{$hash{$t[3]}}, $gene);
}
close(IN);

foreach $id (sort keys %hash){
	@vars = @{$hash{$id}};
	@uvars = keys %{{ map{$_ => 1} @vars}};
	$v = join(",", sort @uvars);
	print "$id\t$v\n";
}