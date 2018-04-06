#!/usr/bin/perl
%distance = ();
%info = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	# $id = $t[0].":".$t[3];
	# if(not exists $info{$id} or abs($t[5]) < $distance{$id}){
	# 	$info{$id} = $_;
	# 	$distance{$id} = abs($t[5]);
	# }
	$id = $t[3].":".$t[9];
	if(not exists $info{$id} or abs($t[12]) < $distance{$id}){
		$info{$id} = $_;
		$distance{$id} = abs($t[12]);
	}
}
close(IN);

foreach $id (keys %info){
	print "$info{$id}\n";
}