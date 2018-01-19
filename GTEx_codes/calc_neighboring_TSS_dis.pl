#!/usr/bin/perl
# %tss = ();
# %tssstrand = ();
%varpromoter = ();
# %varpromoterstrand = ();

open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	if($t[10] eq "promoter"){
		# $varpromoterstrand{$t[3]} = $t[9];
		if($t[9] eq "+"){$varpromoter{$t[3]} = $t[1];}
		else{$varpromoter{$t[3]} = $t[2];}
	}
	# else{
	# 	# $tssstrand{$t[3]} = $t[9];
	# 	if($t[9] eq "+"){$tss{$t[3]} = $t[1];}
	# 	else{$tss{$t[3]} = $t[2];}
	# }
}
close(IN);

print "TSS.distance\n";
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	if($t[9] eq "+"){$tss = $t[1];}
	else{$tss = $t[2];}
	$dis = $tss - $varpromoter{$t[3]};
	print "$dis\t$_\n";
}
close(IN);
