#!/usr/bin/perl
use strict;

my $gz   = $ARGV[0];
my $type = $ARGV[1];

die "No file passed [$gz] ... try again ... "        unless ($gz);
die "Invalid file [$gz] ... try again ... "          unless ($gz);
die 'No type used (TSS|TES|Gene) ... try again ... ' unless ($type);
die 'Invalid type (TES|TSS|Gene) ... try again ... ' unless ($type =~ /(TSS|TES|Gene)/);

open my $PIPE, '-|', "zcat $gz";

while (<$PIPE>){
	chomp;
	my ($bin, $name, $chr, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $nExons, $exStarts, $exEnds, $xa, $xb, $xc, $xd, $xe) = (split "\t") ;

	if ($type eq 'TSS'){
		print join("\t", $chr, $txStart - 1, $txStart, $name, $name, $strand)."\n" if ($strand eq "+");
		print join("\t", $chr, $txEnd   - 1, $txEnd  , $name, $name, $strand)."\n" if ($strand eq "-");
	}

	if ($type eq 'TES'){
		print join("\t", $chr, $txEnd   - 1, $txEnd  , $name, $name, $strand)."\n" if ($strand eq "+");
		print join("\t", $chr, $txStart - 1, $txStart, $name, $name, $strand)."\n" if ($strand eq "-");
	}

	if ($type eq 'Gene'){
		print join("\t", $chr, $txStart,   $txEnd,  $name, $name, $strand)."\n";
	}
}
