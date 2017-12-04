#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
<IN>;
while(<IN>){
	chomp;
	@t = split(/\s/);
	@tt = split(/;/, $t[4]);
	%gene = ();
	foreach $p(@tt){
		@ttt = split(/\|/, $p);
		$gene{$ttt[0]} = 1;
	}
	foreach $g(keys %gene){
		$id = join(" ", $t[1], $t[2], $t[3], $g);
		$hash{$id} = $t[0];
	}
}
close(IN);
foreach $k(keys %hash){
	print "$hash{$k}\s$k\n";
}