use strict; 

use Getopt::Long;
use File::Basename; 
use File::Find;
use File::Which;

GetOptions ('t=s'     		=> \(my $chipBG),
	    	'c=s'      		=> \(my $inputBG),
	    	'sT=s'  		=> \(my $chipBGstat),
	    	'sC=s'  		=> \(my $inputBGstat),
	    	'out=s'    		=> \(my $BGout),
	    	'faidx=s' 		=> \(my $fai),
	    	'tmpDir=s' 		=> \(my $tmpDir = '.'),
	    	's=s'      		=> \(my $scoretype = 'KmetStat'));

die unless (-e $chipBG);
die unless (-e $inputBG);
die if     (-e $BGout);
die unless (-e $fai);
die ('bedGraphToBigWig from UCSC not available!!') unless (which('bedGraphToBigWig'));

my $BGKmet = $BGout; $BGKmet =~ s/(bg|bedgraph)$//; $BGKmet = $BGKmet.'KMETChromosomes.bedgraph';

die unless (-e $chipBGstat);
die unless (-e $inputBGstat);

my (%chipVal,%inputVal);

getKmetStatVals($chipBGstat,\%chipVal);
getKmetStatVals($inputBGstat,\%inputVal);

my $H4me3A_IPe = $chipVal{'KmetStat_H3K4me3_A'} / $inputVal{'KmetStat_H3K4me3_A'};
my $H4me3B_IPe = $chipVal{'KmetStat_H3K4me3_B'} / $inputVal{'KmetStat_H3K4me3_B'};
my $H4me2A_IPe = $chipVal{'KmetStat_H3K4me2_A'} / $inputVal{'KmetStat_H3K4me2_A'};
my $H4me2B_IPe = $chipVal{'KmetStat_H3K4me2_B'} / $inputVal{'KmetStat_H3K4me2_B'};

my $IPe_H3K4me3 = ($chipVal{'KmetStat_H3K4me3_A'} + $chipVal{'KmetStat_H3K4me3_B'}) / ($inputVal{'KmetStat_H3K4me3_A'} + $inputVal{'KmetStat_H3K4me3_A'});
my $nFI         = ($inputVal{'KmetStat_H3K4me3_A'} + $inputVal{'KmetStat_H3K4me3_A'})/($inputVal{'tot'});
my $nFC         = ($chipVal{'KmetStat_H3K4me3_A'} + $chipVal{'KmetStat_H3K4me3_A'})/($chipVal{'tot'});
my $nF          = $nFI/$nFC;

my $tmpBGInit   = $tmpDir.'/init.bg';
my $tmpBGCorr   = $tmpDir.'/corrected.bg';

my $tmpBGKM     = $BGKmet;

open BGOUT,  '>', $tmpBGInit;
open BGKMET, '>', $tmpBGKM;

my %HMDcs;

open my $BGPIPE, '-|', "paste $chipBG $inputBG";
#open BGChIP, $chipBG;
#open BGInput, $inputBG;

#while (defined(my $nChIP=<BGChIP>) && defined(my $nInput=<BGInput>)){
while (<$BGPIPE>){
	chomp ;
	my ($cs,$from,$to,$cVal,$mt1,$mt2,$mt3,$iVal) = split(/\t/,$_);

	my $HMD;

	if ($iVal == 0 || $cVal == 0){
		$HMD=0;
	}else{
		$HMD = sprintf("%4.3f",100 * (($cVal/$iVal)/$IPe_H3K4me3));

		$HMDcs{$cs}->{val} += $HMD;
		$HMDcs{$cs}->{cnt} ++;
	}

	if ($cs =~ /Kmet/){
		print BGKMET join("\t",$cs,$from,$to,$HMD)."\n";
	}else{
		print BGOUT  join("\t",$cs,$from,$to,$HMD)."\n" ;
	}
}

close BGOUT; 
close BGKMET;

normByCSmean($tmpBGInit,\%HMDcs,$BGout);

#########################################################################
sub getKmetStatVals{
	my ($b,$ret) = @_;
	open IN, $b;
	while (<IN>){
		chomp;
		my @F = split(/\t/,$_);
		$$ret{$F[0]}      = $F[1];
		$$ret{'totKmet'} += $F[1] if ($_ =~ /Kmet/);
		$$ret{'tot'}      = $F[2];
	}
	close IN;
}

#########################################################################
sub normByCSmean{
	my ($inBG,$csMean,$outBG) = @_;
	open INBG, $inBG;
	open OUTBG, '>' ,$outBG;
	while (<INBG>){
		chomp;
		my @F = split(/\t/,$_);
		if ($F[0] =~ /^chr\S+$/){
			my $normHMD = $$csMean{$F[0]}->{val}/$$csMean{$F[0]}->{cnt};
			print OUTBG join("\t",@F[0..2],($F[3]?($F[3]-$normHMD):0)."\n");
		}else{
			print OUTBG join("\t",@F[0..3]."\n");
		}
	}
	close INBG;close OUTBG;
}





























































