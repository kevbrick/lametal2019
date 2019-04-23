#!/usr/bin/perl
use strict;
use Getopt::Long;
use Math::Round;

GetOptions ('fimo=s'  => \(my $fimo),
			'w=i'     => \(my $nWin = 250),
			'out=s'   => \(my $out = 'HS.singleMotif.bed'));

my @testPVals = reverse(1e-3,5e-3,1e-4,2e-4,4e-4,6e-4,8e-4,1e-5,2e-5,4e-5,6e-5,8e-5,1e-6,1e-7,1e-8);

my (%hsCnt);

## Parse FIMO file and count motifs at hotspots for each p-value threshold
open(IN, $fimo);

while (<IN>){
	chomp;
	next if ($_ =~ /seque/);
	my ($mName,$mAlt,$hsName,$from,$to,$strand,$score,$P,$Q,$seqMatch) = split(/\t/,$_);
	my ($hsCS,$hsFrom,$hsTo) 					 	  			 = split(/_/,$hsName);
	for my $testP(@testPVals){
		if ($P <= $testP) {
			$hsCnt{$hsName}->{$testP}->{N}++;
			$hsCnt{$hsName}->{$testP}->{fimo} = $_;
		}
	}	
}
close IN;

## Count number of unique hotspots with a single motif for each P
my (%UniqueCountsByP,$totalHS);
for my $hsName(sort keys(%hsCnt)){
	$totalHS++;
	for my $testP(@testPVals){
		next unless ($hsCnt{$hsName}->{$testP}->{N});
		$UniqueCountsByP{$testP}++ if ($hsCnt{$hsName}->{$testP}->{N} == 1);
	}
}

## Get best P-value
my $nMax = -1;
my $bestP;

for my $testP(sort {$a <=> $b} keys(%UniqueCountsByP)){
	if ($UniqueCountsByP{$testP} > $nMax){
		$bestP = $testP;
		$nMax  = $UniqueCountsByP{$testP};
	}
}

## Write BED of HS with one motif
## recenter to motif
open(TMP, '>', 'uniqueHS.tmp');

for my $hsName(sort keys(%hsCnt)){
	if ($hsCnt{$hsName}->{$bestP}->{N} == 1){
		
		my $data =  $hsCnt{$hsName}->{$bestP}->{fimo};

		my ($mName,$mAlt,$hsName,$from,$to,$strand,$score,$P,$Q,$seqMatch) = split(/\t/,$data);
		my ($hsCS,$hsFrom,$hsTo) 					 	  			       = split(/_/,$hsName);

		if ($strand eq "+"){
			$from = $hsFrom + ($from-$nWin-1);
			$to   = $hsFrom + $nWin;
		}else{
			$from = $hsFrom + $to - $nWin;
			$to   = $hsFrom + $to + $nWin + 1;
		}

		print TMP join("\t",$hsCS,$from,$to,$score,$P,$strand)."\n";
	}
}
close TMP;

system('sort -k1,1 -k2n,2n -k3n,3n uniqueHS.tmp >'.$out)

