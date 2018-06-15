#!/usr/bin/perl
### GTRD cell line count from bedtools intersect
%maxcount = ();
%genes = ();
%tfs = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$tfs{$t[3]} = 1;
	$genes{$t[9]} = 1;
	$id = $t[3].":".$t[9];
	if(not exists $maxcount{$id}){
		$maxcount{$id} = $t[4];
	}
	elsif($t[4] > $maxcount{$id}){
		$maxcount{$id} = $t[4];
	}
}
close(IN);

print "Gene";
foreach $tf(sort keys %tfs){
	print "\t$tf";
}
print "\n";

foreach $gene(sort keys %genes){
	print "$gene";
	foreach $tf(sort keys %tfs){
		$id = $tf.":".$gene;
		if(not exists $maxcount{$id}){
			print "\t0";
		}
		else{
			print "\t$maxcount{$id}";
		}
	}
	print "\n";
}
