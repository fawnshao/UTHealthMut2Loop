#!/usr/bin/perl
# initialization
my %motif;
# read in bedtools results
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	my @t = split(/\t/);
	$mkey = join("\t", $t[0], $t[1], $t[2], $t[4], $t[5]);
	unless (exists $motif{$mkey}) {
		$motif{$mkey} = [];
	}
	push (@{$motif{$mkey}}, $t[9]);
}
close(IN);

# count
foreach my $id (keys %motif){
	# my @names = @{$motif{$id}};
	my @tfs = sort keys { map { $_ => 1 } @{$motif{$id}} };
	my $tflist = join(", ", @tfs);
	my $size = scalar @tfs;
	print "$id\t$size\t$tflist\n";
}

