#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	@tt = split(/;/, $t[1]);
	foreach $g(@tt){
		@ttt = split(/\|/, $g);
		unless (exists $hash{$t[0]}) {
			$hash{$t[0]} = [];
		}
		push (@{$hash{$t[0]}}, $ttt[0]);
	}
}
close(IN);
foreach my $k(keys %hash){
	my @genes = sort keys { map { $_ => 1 } @{$hash{$k}} };
	my $gs = join(",", @genes);
	print "$k\t$gs\n";
}