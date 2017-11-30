#!/bin/perl
%exp = ();
%mut = ();
%expid = ();
%mutid = ();
%expgene = ();
%mutpos = ();
%both = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[0].":".$t[1];
	$exp{$id} = $t[2];
	$expid{$t[0]} = 1;
	$expgene{$t[1]} = 1;
}
close(IN);

open(IN, $ARGV[1]) or die "can not open $ARGV[1]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[3].":".$t[0].":".$t[1].":".$t[2];
	$mut{$id} = 1;
	$mutid{$t[3]} = 1;
	$id = $t[0].":".$t[1].":".$t[2];
	$mutpos{$id} = 1;
}
close(IN);

foreach $es (sort keys %expid){
	if (exists $mutid{$es}) {
		$both{$es} = 1;
	}
}

$expout = $ARGV[2].".bothWGS.exp.tsv";
open(OUT, ">$expout");
print OUT "Gene";
foreach $sample(sort keys %both){
	print OUT "\t$sample";
}
print OUT "\n";

foreach $g(sort keys %expgene){
	print OUT $g;
	foreach $sample(sort keys %both){
		$id = $sample.":".$g;
		if(not exists $exp{$id}){
			$exp{$id} = 0
		}
		print OUT "\t$exp{$id}";
	}
	print OUT "\n";
}
close(OUT);

$mutout = $ARGV[2].".bothWGS.mut.tsv";
open(OUT, ">$mutout");
print OUT "Mut";
foreach $sample(sort keys %both){
	print OUT "\t$sample";
}
print OUT "\n";

foreach $g(sort keys %mutpos){
	print OUT $g;
	foreach $sample(sort keys %both){
		$id = $sample.":".$g;
		if(not exists $mut{$id}){
			$mut{$id} = 0
		}
		print OUT "\t$mut{$id}";
	}
	print OUT "\n";
}
close(OUT);
