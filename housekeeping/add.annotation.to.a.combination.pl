#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]);
while (<IN>) {
	chomp;
	@t = split(/\t/);
	$hash{$t[0]} = $t[1];
}
close(IN);

open(IN, $ARGV[1]);
while (<IN>) {
	chomp;
	@t = split(/\t/);
	$a = "";
	$b = "";
	if($t[2] ne "na"){
		@tt = split(/,/, $t[2]);
		foreach $x(@tt){
			if(not exists $hash{$x}){
				$hash{$x} = "null";
			}
			$a = $hash{$x}.",";
		}
		$a =~ s/,$//;
	}
	if($t[3] ne "na"){
		@tt = split(/,/, $t[3]);
		foreach $x(@tt){
			if(not exists $hash{$x}){
				$hash{$x} = "null";
			}
			$b = $hash{$x}.",";
		}
		$b =~ s/,$//;
	}
	print "$_\t$a\t$b\n";
}
close(IN);
