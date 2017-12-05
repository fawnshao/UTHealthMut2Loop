#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
$line = <IN>;
chomp $line;
print "$line\n";
while(<IN>){
	chomp;
	@t = split(/\s/);
	@tt = split(/,/, $t[4]);
	%gene = ();
	foreach $p(@tt){
		@ttt = split(/;/, $p);
		foreach $pp(@ttt){
			@tttt = split(/\|/, $pp);
			$gene{$tttt[0]} = 1;
		}
	}
	$mutp="";
	foreach $g(keys %gene){
		$mutp .= $g.",";
	}
	$mutp =~ s/,$//;
	$id = join(" ", $t[1], $t[2], $t[3], $mutp);
	$hash{$id} = $t[0];
}
close(IN);
foreach $k(keys %hash){
	print "$hash{$k} $k\n";
}