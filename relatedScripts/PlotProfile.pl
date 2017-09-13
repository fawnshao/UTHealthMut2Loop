#!/usr/bin/perl
$rfile = "/home1/04935/shaojf/myTools/UTHealthMut2Loop/relatedScripts/PlotProfile.R";
if(@ARGV < 3) {
	die "Usage: perl $0 <input bedtools intersect file> <y for directionality> <output>\n";
}

open(OUT, ">$ARGV[2]") or die "can not open $ARGV[2]\n";
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	$linesize = @t;
	if($t[$linesize - 1 ] > 0){
		$mid = $t[1] + 2000;
		$petmid = ($t[$linesize - 5] + $t[$linesize - 6]) / 2;
		@a = split(/,/, $t[$linesize - 4]);
		$a[0] =~ s/\.\./:/g;
		$a[0] =~ s/-/:/g;
		@b = split(/:/, $a[0]);
		$leftmid = ($b[1] + $b[2]) / 2;
		$rightmid = ($b[4] + $b[5]) / 2;
		if(abs($leftmid - $petmid) < 5){
			$left = $b[1];
			$right = $b[2];
		}
		else{
			$left = $b[4];
			$right = $b[5];
		}
		if($ARGV[1] eq 'y' && $t[5] eq '-'){
			$pos2 = $mid - $left;
			$pos1 = $mid - $right;
		}
		else{
			$pos1 = $left - $mid;
			$pos2 = $right - $mid;
		}
		print OUT "$pos1\t$pos2\t$a[1]\n";
	}
}
close(IN);
close(OUT);

$cmd = "module load Rstats; Rscript $rfile $ARGV[2]";
system($cmd);
