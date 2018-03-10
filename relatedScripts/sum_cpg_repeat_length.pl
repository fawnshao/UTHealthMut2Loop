#!/usr/bin/perl
### CpG island or simple repeat length sum from bedtools intersect
%sum = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$sum{$t[7]} += $t[10];
}
close(IN);

print "Gene\t$ARGV[0]\n";

foreach $gene(sort keys %sum){
	print "$gene\t$sum{$gene}\n";
}
