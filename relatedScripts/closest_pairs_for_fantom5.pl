#!/usr/bin/perl
# bedtools closest file: fantom5.v19.promoter.within1k.*txt
# grep -v "@chr" fantom5.v19.promoter.within1k.div.txt | grep -v ENST | grep -v ",p" | wc -l
# cage_peaks annotation file: hg19.cage_peak_phase1and2combined_anno.txt
# entrenz gene annotation: Homo_sapiens.gene_info.annotation.txt 
# exclude cage peaks like p1@ZNF10,p2@ZNF268?
# $cage2entrezgene = ();
$entrezgene2type = ();

# Homo_sapiens.gene_info.annotation.txt
open(IN, $ARGV[1]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$entrezgene2type{$t[1]} = $t[2];
}
close(IN);

## bedtools result
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	@tt = split(/@/, $t[3]);
	@ttt = split(/@/, $t[9]);
	if(exists $entrezgene2type{$tt[1]} && exists $entrezgene2type{$ttt[1]}){
		print "$t[3]\t$t[9]\t$t[12]\t$entrezgene2type{$tt[1]}\t$entrezgene2type{$ttt[1]}\n";
	}
}
close(IN);

