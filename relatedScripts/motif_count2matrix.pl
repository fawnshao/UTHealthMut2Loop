#!/usr/bin/perl
if(@ARGV < 2){
	die "usage: perl $0 <input> <motif column> > <output>\nplease use the bedtools intersect results\n";
}

$motif = $ARGV[1];
%hash = ();
%genes = ();
%motifs = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = join("\t", $t[0], $t[1], $t[2], $t[3]);
	$mykey = $id."&&".$t[$motif];
	unless (exists $hash{$mykey}) {
		$hash{$mykey} = 0;
	}
	$hash{$mykey}++;
	$genes{$id} = 1;
	$motifs{$t[$motif]} = 1;
}
close(IN);

print "chr\tstart\tend\tname";
foreach $m (sort keys %motifs){
	print "\t$m";
}
print "\n";
foreach $g (sort keys %genes){
	print "$g";
	foreach $m (sort keys %motifs){
		$mkey = $g."&&".$m;
		if(not exists $hash{$mkey}){$hash{$mkey} = 0;}
		print "\t$hash{$mkey}";
	}
	print "\n";
}