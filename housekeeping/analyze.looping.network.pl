#!/usr/bin/perl
# input a <node1 node2 interation> file
# output degree, subnetwork size
%degrees = ();
%subnetworks = ();
%nodes = ();
$i = 1;
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	if($t[1] =~ /:/){
		@tt = split(/:/, $t[1]);
		$t[1] = $tt[1];
	}
	$degrees{$t[0]}++;
	$degrees{$t[1]}++;
	if(not exists $subnetworks{$t[0]} && exists $subnetworks{$t[1]}){
		$subnetworks{$t[0]} = $subnetworks{$t[1]};
		$nodes{$subnetworks{$t[1]}} .= "\t".$t[0];
	}
	if(exists $subnetworks{$t[0]} && not exists $subnetworks{$t[1]}){
		$subnetworks{$t[1]} = $subnetworks{$t[0]};
		$nodes{$subnetworks{$t[0]}} .= "\t".$t[1];
	}
	if(not exists $subnetworks{$t[0]} && not exists $subnetworks{$t[1]}){
		$subnetworks{$t[0]} = $i;
		$subnetworks{$t[1]} = $i;
		$nodes{$i} = join("\t", $t[0], $t[1]);
		$i++;
	}
}
close(IN);

open(OUT, ">$ARGV[0].nodes.degrees.txt");
foreach $n(keys %degrees){
	print OUT "$n\t$degrees{$n}\n";
}
close(OUT);

open(OUT, ">$ARGV[0].node2subnetworks.txt");
foreach $s(keys %subnetworks){
	print OUT "$s\tsubnet.$subnetworks{$s}\n";
}
close(OUT);

open(OUT, ">$ARGV[0].subnetworks.txt");
foreach $s(keys %nodes){
	@t = split(/\t/, $nodes{$s});
	$size = scalar(@t);
	print OUT "subnet.$s\t$size\t$nodes{$s}\n";
}
close(OUT);
