#!/usr/bin/perl
%tfs = ();
%genes = ();
%hash = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	@tt = split(/\|/, $t[3]);
	@ttt = split(/;/, $tt[1]);
	$genes{$t[7]} = 1;
	foreach $c(@ttt){
		# print "$tt[0]%$c\t$t[7]\n";
		$tf = $tt[0]."%".$c;
		$tfs{$tf} = 1;
		$id = $tf.":".$t[7];
		$hash{$id} = 1;
	}
}
close(IN);

print "Gene";
foreach $a(sort keys %tfs){
	print "\t$a";
}
print "\n";
foreach $g(keys %genes){
	print "$g";
	foreach $a(sort keys %tfs){
		$id = $a.":".$g;
		print "\t$hash{$id}";
	}
	print "\n";
}
