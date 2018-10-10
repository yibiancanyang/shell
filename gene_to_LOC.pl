#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<>){
	chomp;
	if($_=~ /.*\tID=(gene[0-9]*);.*Name=(LOC[0-9]*);.*$/){
		print OUT "$1\t$2\n"; 
		}
	
}
close IN;
close OUT;
