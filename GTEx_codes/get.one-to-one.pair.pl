#!/usr/bin/perl
# awk '$1==1{print $2}' Cells_EBV-transformed_lymphocytes.v7.eVariants.count.txt > single.esnp.txt
# awk '$1==1{print $2}' Cells_EBV-transformed_lymphocytes.v7.eGene.count.txt > single.egene.txt
%singlesnp = ();
%singlegene = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	$singlesnp{$_} = 1;
}
close(IN);

open(IN, $ARGV[1]);
while(<IN>){
	chomp;
	$singlegene{$_} = 1;
}
close(IN);

open(IN, $ARGV[2]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	if(exists $singlesnp{$t[0]} && exists $singlegene{$t[1]}){
		print "$_\n";
	}
}
close(IN);
