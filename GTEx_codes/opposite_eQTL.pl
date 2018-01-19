#!/usr/bin/perl
%pslope = ();
%slopes = ();
%pflag = ();
%info = ();

open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[3]."\t".$t[4];
	$varid = $t[3];
	$info{$id} = $_;
	$slopes{$id} = $t[7];
	$pflag{$id} = "/";
	if($t[6] > -1000 && $t[6] < 1000){
		unless (exists $pslope{$varid}) {
			$pslope{$varid} = [];
		}
		push (@{$pslope{$varid}}, $t[7]);
		$pflag{$id} = "promoter";
	}
}
close(IN);

print "gene.chr\tgene.start\tgene.end\tvar\tgene.id\tgene.name\tdis_var_tss\tslope\tNULL\tgene.strand\tpromoter.flag\texp.flag\n";
foreach $pair (sort keys %info){
	($varid, $geneid) = split(/\t/, $pair);
	@promoterslope = @{$pslope{$varid}};
	$flag = "";
	$v = 0;
	$ps = 0;
	if(@promoterslope == 1){
		$ps = $promoterslope[0];
	}
	else{
		$poscount = 0;
		$negcount = 0;
		foreach $a(@promoterslope){
			if($a < 0 ){$negcount++;}
			if($a > 0 ){$poscount++;}
		}
		if($negcount == @promoterslope){$ps = -1;}
		elsif($poscount == @promoterslope){$ps = 1;}
		else{$ps = 0;}
	}
	if($ps == 0){
		$flag = "contradict";
	}
	else{
		$v = $slopes{$pair} / $ps;
		if($v < 0 ){
			$flag = "opposite";
		}
		else{
			$flag = "/";
		}
	}
	print "$info{$pair}\t$pflag{$pair}\t$flag\n";
}