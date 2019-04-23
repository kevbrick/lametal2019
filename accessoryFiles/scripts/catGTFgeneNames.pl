#!/usr/bin/perl
use strict; 

while (<>){
	chomp;
	
	## Header
	if ($_ =~ /^\#\#/){
		print $_."\n";
		next;
	}

	my @F = split(/\t/,$_);

	my ($ensGID,$ensTID,$havGID,$havTID,$nameGID,$nameTID,$support,$type) ; 
	
	if ($_ =~ /gene_id\s\"(\S+)\"/){$ensGID=$1}; 
	if ($_ =~ /transcript_id\s\"(\S+)\"/){$ensTID=$1}; 

	if ($_ =~ /havana_gene\s\"(\S+)\"/){$havGID=$1}; 
	if ($_ =~ /havana_transcript\s\"(\S+)\"/){$havTID=$1}; 

	if ($_ =~ /gene_name\s\"(\S+)\"/){$nameGID=$1}; 
	if ($_ =~ /transcript_name\s\"(\S+)\"/){$nameTID=$1}; 

	if ($_ =~ /transcript_support_level\s\"(\S+)\"/){$support=$1}; 
	if ($_ =~ /transcript_type\s\"(\S+)\"/){$type=$1}; 

	my $myID = join("|",$ensTID,$ensGID,$havTID,$havGID,$nameTID,$nameGID,$support,$type); 

	$_ =~ s/transcript_id\s\"(\S+?)\"/transcript_id "$myID"/; 

	print $_ ."\n" ;
}
