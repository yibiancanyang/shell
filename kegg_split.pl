#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<IN>){
	chomp;
	my @array=split(/\t/);
	if($array[1]=~ /.*(KO\:K[0-9]*)/){
		print OUT "$array[0]\t$1\n"; 
		}
	
}
close IN;
close OUT;

