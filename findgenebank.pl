#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<>){
	chomp;
	if($_=~ /.*\t(ID=rna[0-9]*);.*(Genbank:[A-Z0-9_]*\.1);.*/){
		print OUT "$1\t$2\n"; 
		}
	
}
close IN;
close OUT;

