use strict; 

use Getopt::Long;
use File::Basename; 
use File::Find;
use File::Which;

GetOptions ('t=s'     		=> \(my $chipBG),
	    	'c=s'      		=> \(my $inputBG),
	    	'ncis=s'   		=> \(my $ncis),
	    	'out=s'    		=> \(my $BGout),
	    	'tmpDir=s' 		=> \(my $tmpDir = '.'));

if ($BGout !~ /bedgraph$/){
	$BGout .= '.bedgraph';
}

die unless (-e $chipBG);
die unless (-e $inputBG);
die if     (-e $BGout);

my $tmpBG       = $tmpDir.'/bg_'.$BGout.'.tmp';

open BGOUT,  '>', $BGout;

open BGChIP, $chipBG;
open BGInput, $inputBG;

while (defined(my $nChIP=<BGChIP>) && defined(my $nInput=<BGInput>)){
	chomp $nChIP; chomp $nInput;
	my ($cCS,$cFrom,$cTo,$cVal) = split(/\t/,$nChIP);
	my ($iCS,$iFrom,$iTo,$iVal) = split(/\t/,$nInput);

	my $score = sprintf("%4.3f",($cVal - $iVal*$ncis));
	
	print BGOUT  join("\t",$cCS,$cFrom,$cTo,$score)."\n" ;
}

close BGOUT; 