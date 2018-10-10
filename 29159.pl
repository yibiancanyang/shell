#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<IN>){
	chomp;
	if($_=~ /.*\|([A-Z]+_[0-9\.]*)\|.*/){
		print OUT ">$1\n"; 
		}
		else {print OUT "$_\n";}
	
}
close IN;
close OUT;
