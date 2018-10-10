#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<>){
	chomp;
	if($_=~ /.*\tgene_id \"([A-Za-z0-9\.]+)\";.*; gene_name \"(\w+)\";.*$/){
		print OUT "$1\t$2\n"; 
		}
	
}
close IN;
close OUT;
