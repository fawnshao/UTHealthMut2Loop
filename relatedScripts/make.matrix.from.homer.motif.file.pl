#!/usr/bin/perl
### annotatePeaks.pl
# after cut
open(IN, $ARGV[0]);
<IN>;
chomp;
print "$_\n";
while(IN){
	chomp;
	@t = split(/\t/);
	print "$t[0]";
	for($i = 1; $i < @t; $i++){
		@tt = split(/),/, $t[$i]);
		$a = length(@tt);
		print "\t$a";
	}
	print "\n";
}
close(IN);