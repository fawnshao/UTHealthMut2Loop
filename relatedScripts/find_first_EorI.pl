#!/usr/bin/perl
%hash = ();
%pos = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	@tt = split(/:/, $t[3]);
	if(not exists $hash{$tt[0]}){
		$hash{$tt[0]} = $_;
		$pos{$tt[0]} = $t[1];
	}
	else{
		if($t[5] eq "+" && $t[1] < $pos{$tt[0]}){
			$hash{$tt[0]} = $_;
			$pos{$tt[0]} = $t[1];
		}
		elsif($t[5] eq "-" && $t[1] > $pos{$tt[0]}){
			$hash{$tt[0]} = $_;
			$pos{$tt[0]} = $t[1];
		}
	}
}
close(IN);
foreach $g(keys %hash){
	print "$hash{$g}\n";
}