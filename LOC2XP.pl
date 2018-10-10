#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<IN>){
	chomp;
	if($_=~ /Name=(XP_[0-9\.]*);.*gbkey=CDS;gene=([A-Za-z0-9_-]*);/){
		print OUT "$1\t$2\n"; 
		}
	
}
close IN;
close OUT;
