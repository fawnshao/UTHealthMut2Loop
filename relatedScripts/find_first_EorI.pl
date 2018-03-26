#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	@tt = split(/:/, $t[3]);
	if($t[5] eq "+" && not exists $hash{$tt[0]}){
		$hash{$tt[0]} = $_;
	}
	elsif($t[5] eq "-"){
		$hash{$tt[0]} = $_;
	}
}
close(IN);
foreach $g(keys %hash){
	print "$hash{$g}\n";
}