#!/usr/bin/perl
### annotatePeaks.pl
# after cut
open(IN, $ARGV[0]);
$line = <IN>;
chomp $line;
@t = split(/\t/,$line);
print "$t[0]";
for($i = 1; $i < @t; $i++){
	print "\t$t[$i]";
}
print "\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	print "$t[0]";
	print STDERR "@t\n";
	for($i = 1; $i < @t; $i++){
		if($t[$i] eq ''){
			$a = 0;
		}
		else{
			@tt = split(/\),/, $t[$i]);
			$a = @tt;
		}
		print "\t$a";
	}
	print "\n";
}
close(IN);