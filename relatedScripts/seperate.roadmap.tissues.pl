#!/usr/bin/perl
%tfs = ();
%genes = ();
%hash = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	@tt = split(/_/, $t[3]);
	$tf = $tt[0]."_".$tt[1];
	$tfs{$tf} = 1;
	$id = $tf.":".$t[7];
	$hash{$id} = 1;
	$genes{$t[7]} = 1;
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
		if(exists $hash{$id}){
			print "\t$hash{$id}";
		}
		else{
			print "\t0";
		}
	}
	print "\n";
}
