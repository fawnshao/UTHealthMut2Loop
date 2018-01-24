#!/usr/bin/perl
$flag = ();
$keynum = $ARGV[2];
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	$flag{$_} = 1;
}
close(IN);

open(IN, $ARGV[1]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	if(exists $flag{$t[$keynum]}){
		print "$_\n";
	}
}
close(IN);