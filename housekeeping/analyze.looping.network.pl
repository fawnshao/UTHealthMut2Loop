#!/usr/bin/perl
# input a <node1 node2 interaction> file
# output degree, subnetwork size
%degrees = ();
%subnetworks = ();
%nodes = ();
####### add a relationship between interaction and subnetwork
%interaction = ();
#######
$i = 0;
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
<IN>;
while(<IN>){
	chomp;
	@t = split(/\t/);
	# if($t[1] =~ /:/){
	# 	@tt = split(/:/, $t[1]);
	# 	$t[1] = $tt[1];
	# }
	$degrees{$t[0]}++;
	$degrees{$t[1]}++;
	if(! exists $subnetworks{$t[0]} && ! exists $subnetworks{$t[1]}){
		$i++;
		$subnetworks{$t[0]} = $i;
		$subnetworks{$t[1]} = $i;
		$nodes{$i} = $t[0]."\t".$t[1];
	}
	elsif(! exists $subnetworks{$t[0]} && exists $subnetworks{$t[1]}){
		$subnetworks{$t[0]} = $subnetworks{$t[1]};
		$nodes{$subnetworks{$t[1]}} .= "\t".$t[0];
	}
	elsif(exists $subnetworks{$t[0]} && ! exists $subnetworks{$t[1]}){
		$subnetworks{$t[1]} = $subnetworks{$t[0]};
		$nodes{$subnetworks{$t[0]}} .= "\t".$t[1];
	}
	##### for the above processing, we would miss a situation which is:
	##### both of the nodes have already been assigned to different subnetworks
	##### that we need to merge them
	##### for v2 modification, we will try to take this into consideration!!!!
	elsif(exists $subnetworks{$t[0]} && exists $subnetworks{$t[1]} && $subnetworks{$t[0]} != $subnetworks{$t[1]}){
		$nodes{$subnetworks{$t[0]}} .= "\t".$nodes{$subnetworks{$t[1]}};
		@x = split(/\t/,$nodes{$subnetworks{$t[1]}});
		delete $nodes{$subnetworks{$t[1]}};
		foreach $xx(@x){
			$subnetworks{$xx} = $subnetworks{$t[0]};
		}
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

open(OUT, ">$ARGV[0].interaction2subnetworks.txt");
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
<IN>;
while(<IN>){
	chomp;
	@t = split(/\t/);
	$intnodes = $t[0]."\t".$t[1];
	$interaction{$intnodes} = "subnet".$subnetworks{$t[0]};
	print OUT "$intnodes\t$interaction{$intnodes}\n";
}
close(IN);
close(OUT);
