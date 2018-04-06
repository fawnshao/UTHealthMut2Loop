#!/usr/bin/perl
%distance = ();
%info = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[0].":".$t[3];
	if(not exists $info{$id} or abs($t[5]) < $distance{$id}){
		$info{$id} = $_;
		$distance{$id} = abs($t[5]);
	}
}
close(IN);

foreach $id (keys %info){
	print "info{$x}\n";
}