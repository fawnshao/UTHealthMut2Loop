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
	if(not exists)
}
close(IN);