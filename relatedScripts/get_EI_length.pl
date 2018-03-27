#!/usr/bin/perl
%gene = ();
%intron = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$gene{$t[3]} = $t[2] - $t[1] + 1;
}
close(IN);

open(IN, $ARGV[1]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	if(not exists $intron{$t[3]}){
		$intron{$t[3]} = 0;
	}
	$intron{$t[3]} += $t[2] - $t[1] + 1;
}
close(IN);

foreach $g(keys %gene){
	if(not exists $intron{$g}){
		$intron{$g} = 0;
	}
	print STDERR "$g\t$intron{$g}\n";
	$x = $gene{$g} - $intron{$g};
	print "$g\t$x\n";
}