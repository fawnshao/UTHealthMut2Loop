#!/usr/bin/perl
%motifs = ();
%count = ();
open(IN, $ARGV[0]);
$line = <IN>;
chomp $line;
@t = split(/\t/, $line);
$max = scalar(@t);
# print STERR "$max\n";
for($i = 5; $i < $max; $i++){
	@tt = split(/\//, $t[$i]);
	$motifs{$i} = $tt[0];
	# print STERR "$i\t$motifs{$i}\n";
}
while(<IN>){
	chomp;
	@t = split(/\t/);
	for($i = 5; $i < $max; $i++){
		if($t[$i] =~ /,/){
			@tt = split(/;/, $t[$i]);
			foreach $a(@tt){
				@ttt = split(/\(/, $a);
				$b = $ttt[0];
				# print STERR "$b\n";
				if($t[4] eq '-'){
					$dis = (1000 - $b); #/ 100;
				}
				else{
					$dis = ($b - 1000); #/ 100;
				}
				$bin = sprintf("%0d", $dis);
				$id = $motifs{$i}."\t".$bin;
				# print STERR "$id\n";
				$count{$id}++;
			}
		}
	}
}
close(IN);

foreach $x (keys %count){
	print "$x\t$count{$x}\n";
}