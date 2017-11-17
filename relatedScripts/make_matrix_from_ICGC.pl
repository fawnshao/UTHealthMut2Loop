#!/usr/bin/perl
%samples = ();
%genes = ();
%exprs = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
<IN>;
while(<IN>){
	chomp;
	@t = split(/\t/);
	$samples{$t[0]} = 1;
	$genes{$t[1]} = 1;
	$id = join("#", $t[0], $t[1]);
	# $exprs{$id} = $t[2];
	# some samples have replications, use the average
	unless (exists $exprs{$id}) {
		$exprs{$id} = [];
	}
	push (@{$exprs{$id}}, $t[2]);
}
close(IN);

foreach $g(keys %genes){
	print "$g";
	foreach $s(keys %samples){
		$id = join("#", $s, $g);
		$v = 0;
		if(exists $exprs{$id}){
			$v = average(@{$exprs{$id}});
		}
		print "\t$v";
	}
	print "\n";
}

sub average{
	my $total;
	$total += $_ for @_;
	# return sprintf '%.2f', $total / @_;
	return $total / @_;
}