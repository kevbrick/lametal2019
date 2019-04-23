#!/usr/bin/perl
use strict;
use Getopt::Long;
use Math::Round;

GetOptions ('e=s' => \(my $expData),
			't=s' => \(my $tssBed));

my %tData;
open INTSS, $tssBed;
while (<INTSS>){
	chomp;
	my @F = split(/\t/,$_);
	$tData{$F[3]}->{cs} = $F[0];
	$tData{$F[3]}->{from} = $F[1];
	$tData{$F[3]}->{to} = $F[2];
}
close INTSS;

print join("\t",'cs','from','to','names','LE','ZY','PA','DI')."\n";

open INEXP, $expData;
while (<INEXP>){
	chomp;
	my ($names,$LE,$ZY,$PA,$DI) = split(/\t/,$_);

	my @N = split(/\|/,$names);

	my $isoform = $N[4];

	die unless($tData{$isoform}->{to});
	print join("\t",$tData{$isoform}->{cs},$tData{$isoform}->{from},$tData{$isoform}->{to},$_)."\n";
}
close INEXP;
