#!/usr/bin/perl
if(@ARGV < 2){
	die "usage: perl $0 <patient-based tsv input> <promoter motif count> > <output>\n\n";
}

%hash = ();
%genes = ();
%motifs = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = join(":", $t[0], $t[1], $t[2], $t[3]);
	$mykey = $id."&&".$t[$motif];
	unless (exists $hash{$mykey}) {
		$hash{$mykey} = 0;
	}
	$hash{$mykey}++;
	$genes{$id} = 1;
	$motifs{$t[$motif]} = 1;
}
close(IN);

print "patient\tdisease\tTAD\tGene\tmut.mean\tctr.mean\tctr.sd\tfc\toutlier.flag\tzscore.1\tzscore.2\tmut.exp\tctr.exp\tmut.flag\talt.flag\tpromoter.motif.count\tpromoter.motifs\n";
foreach $g (sort keys %genes){
	print "$g";
	foreach $m (sort keys %motifs){
		$mkey = $g."&&".$m;
		if(not exists $hash{$mkey}){$hash{$mkey} = 0;}
		print "\t$hash{$mkey}";
	}
	print "\n";
}