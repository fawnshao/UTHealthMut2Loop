#!/usr/bin/perl
%mut = ();
%ref = ();
%all = ();

$refin = $ARGV[0].".ref.fimo.txt";
$mutin = $ARGV[0].".mut.fimo.txt";
# $gain = $ARGV[0].".LoopBroken.motif.gain";
# $lose = $ARGV[0].".LoopBroken.motif.lose";
$out = $ARGV[0].".motif.cmp";

open(IN, $refin);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[0]."\t".$t[1];
	$ref{$id} = 1;
	$all{$id} = 1;
}
close(IN);

open(IN, $mutin);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[0]."\t".$t[1];
	$mut{$id} = 1;
	$all{$id} = 1;
}
close(IN);

%gain = ();
%lose = ();
%var = ();
foreach $pair (sort keys %all){
	($motif, $varid) = split(/\t/, $pair);
	if($ref{$pair} == 1 && $mut{$pair} != 1){
		$lose{$varid} .= $motif . ","; 
		$var{$varid} = 1;
	}
	if($ref{$pair} != 1 && $mut{$pair} == 1){
		$gain{$varid} .= $motif . ","; 
		$var{$varid} = 1;
	}
}

open(OUT, ">$out");
print OUT "var\tgain\tlose\n";
foreach $v (sort keys %var){
	if($v ne ""){
		if(not exists $gain{$v}){$gain{$v} = "/";}
		if(not exists $lose{$v}){$lose{$v} = "/";}
		$gain{$v} =~ s/,$//;
		$lose{$v} =~ s/,$//;
		print OUT "$v\t$gain{$v}\t$lose{$v}\n";
	}
}
close(OUT);

