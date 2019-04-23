#!/usr/bin/perl
use strict; 

system('awk \'{close(f);f=$1}{print > f".bed"}\' '.$ARGV[0]);

my $catCmd = 'cat ';
my $rmCmd  = 'rm ';

open CS, $ARGV[1];

while (<CS>){
	chomp; 
	my @F = split(/\t/,$_);
	if (-s $F[0].'.bed'){
		system("sort -k1,1 -k2n,2n -k3n,3n $F[0]\.bed >$F[0].sorted.bed");
		$catCmd .= " $F[0].sorted.bed";
		$rmCmd  .= " $F[0].sorted.bed";
	}
	$rmCmd .= " $F[0].bed" if (-e "$F[0].bed");
}

close CS; 

system($catCmd);
system($rmCmd);
