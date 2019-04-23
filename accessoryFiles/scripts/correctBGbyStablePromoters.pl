use strict; 

use Getopt::Long;
use File::Which;

GetOptions ('bg=s'    	=> \(my $bedgraph),
			'tss=s'     => \(my $stableTSS),
	    	'out=s'   	=> \(my $BGout),
	    	'tmpDir=s'	=> \(my $tmpDir = '.'));

die unless (-e $bedgraph);
die unless (-e $stableTSS);
die if     (-e $BGout);

my $tmpBG       = $tmpDir.'/bg_'.$BGout.'.tmp';

#my $stableTSS = ($bedgraph =~ /H3K4me3/)?$stableK4:$stableK9;

open my $PIPE, '-|', "intersectBed -a $bedgraph -b $stableTSS |cut -f4";

my ($sum,$tot,$mean);

while (<$PIPE>){
	chomp; 
	$tot++;
	$sum+=$_;
	$mean=$sum/$tot;
}

close $PIPE;

open BGOUT,  '>', $BGout;

open BG, $bedgraph;

while (<BG>){
	chomp;

	my ($cs,$from,$to,$v) = split(/\t/,$_);
	my $score = sprintf("%4.3f",$v/$mean);
	
	print BGOUT  join("\t",$cs,$from,$to,$score)."\n" ;
}

close BGOUT; 