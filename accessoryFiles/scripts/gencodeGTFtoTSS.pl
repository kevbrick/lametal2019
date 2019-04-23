#!/usr/bin/perl

use strict; 

my $gtf = $ARGV[0];

open my $PIPE, '-|', ($gtf =~ /gz$/)?"zcat $gtf":"cat $gtf";

while (<$PIPE>){
	chomp; 
	my ($chr, $source, $feat, $start, $end, $score, $strand, $frame, $meta) = (split "\t") ;
	next unless ($feat eq "transcript");
	
	my $gtype = ($meta =~ /\"protein_coding\"/)?"proteincoding":"other"; 
	
	my ($name,$gname); 
	$name  = $1 if ($meta =~ /transcript_name\s\"(\S+?)\"/);
	$gname = $1 if ($meta =~ /gene_name\s\"(\S+?)\"/);
	
	print join("\t", $chr, $start,   $start +1, $name, $gtype.":".$gname, $strand)."\n" if ($strand eq "+");
	print join("\t", $chr, $end - 1, $end,      $name, $gtype.":".$gname, $strand)."\n" if ($strand eq "-");
}
