#!/usr/bin/perl
### GTRD cell line count from bedtools intersect
%count = ();
%genes = ();
%tfs = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$tfs{$t[0]} = 1;
	$genes{$t[1]} = 1;
	$id = $t[0].":".$t[1];
	$count{$id} = $t[3];
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
		if(not exists $count{$id}){
			print "\t0";
		}
		else{
			print "\t$count{$id}";
		}
	}
	print "\n";
}
