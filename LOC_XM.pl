#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<IN>){
	chomp;
	if($_=~ /^.*\t[a-z_]*RNA\t.*ID=(rna[0-9]*);.*Genbank:([A-Z0-9_\.]*);.*;gbkey=[a-z_]*RNA;gene=([A-Za-z0-9_-]*);.*$/){
		print OUT "$1\t$2\t$3\n"; 
		}
	
}
close IN;
close OUT;
