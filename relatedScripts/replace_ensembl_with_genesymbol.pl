#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]);
print STDERR "reading $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	$hash{$t[0]} = $t[1];
	}
close(IN);
print STDERR "finished $ARGV[0] hash\n";

open(IN,$ARGV[1]);
print STDERR "reading $ARGV[1]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	@id = split(/\./, $t[1]);
	if(not exists $hash{$id[0]}){print "$_\n";}
	else{
		print "$t[0]\t$hash{$id[0]}\t$t[2]\n"
		}
	}
close(IN);
print STDERR "finished $ARGV[1] hash\n";