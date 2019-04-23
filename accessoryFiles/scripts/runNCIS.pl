#!/usr/bin/perl
use strict;
use Getopt::Long;

GetOptions ('t=s' 	  => \(my $treat),
            'c=s'  	  => \(my $input),
            'lib=s'   => \(my $NCISlibPath),
	    	'out=s'	  => \(my $NCISout),
	    	'tmp=s'   => \(my $tmpDir = "."));

my ($treatFile,$inputFile) = ($treat,$input);

my $NCISscript = "doNCIS.R";

open OUT, '>', $NCISscript;

print OUT 'library("NCIS", lib.loc="'.$NCISlibPath.'")'."\n";
print OUT 'library(\'rtracklayer\')'."\n";
print OUT 'library("ShortRead")'."\n";
print OUT 'res <- NCIS("'.$treatFile.'","'.$inputFile.'","BED")'."\n";
print OUT 'write(paste(res$est,res$pi0,res$binsize.est,res$r.seq.depth,sep = "\t"), "'.$NCISout.'", sep = "\t")'."\n";

close OUT;

system('R --silent --vanilla <'.$NCISscript);