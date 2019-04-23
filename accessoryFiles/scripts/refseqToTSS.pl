#!/usr/bin/perl
use strict; 

my $gz = $ARGV[0];

open my $PIPE, '-|', "zcat $gz";

while (<$PIPE>){
	chomp; 
	my ($bin, $name, $chr, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $nExons, $exStarts, $exEnds, $xa, $xb, $xc, $xd, $xe) = (split "\t") ;
	
	print join("\t", $chr, $txStart,   $txStart +1, $name, $name, $strand)."\n" if ($strand eq "+");
	print join("\t", $chr, $txEnd - 1, $txEnd     , $name, $name, $strand)."\n" if ($strand eq "-");
}
