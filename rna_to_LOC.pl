#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<>){
	chomp;
	if($_=~ /.*\tmRNA\t.*ID=(rna[0-9]*);.*;gene=([A-Za-z0-9]*);.*;transcript_id=([A-Za-z0-9_\.]*)/){
		print OUT "$1\t$2\t$3\n"; 
		}
	
}
close IN;
close OUT;