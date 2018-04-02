#!/usr/bin/perl
%motifs = ();
%count = ();
open(IN, $ARGV[0]);
$line = <IN>;
chomp $line;
@t = split(/\t/, $line);
$max = scalar(@t);
for($i = 5; $i < $max; $i++){
	@tt = split(/\//, $t[$i]);
	$motifs{$i} = $tt[0];
}
while(<IN>){
	chomp;
	@t = split(/\t/);
	for($i = 5; $i < $max; $i++){
		if($t[$i] ne ""){
			@tt = split(/,/, $t[$i]);
			foreach $a(@tt){
				@ttt = split(/(/, $a);
				foreach $b(@ttt){
					if($t[4] eq '-'){
						$dis = (2000 - $b) / 100;
					}
					else{
						$dis = $b;
					}
					$bin = sprintf("%0d", $dis);
					$id = $motifs{$i}.":".$bin;
					$count{$id}++;
				}
			}
		}
	}
}
close;

foreach $x (keys %count){
	print "$x\t$count{$x}\n";
}