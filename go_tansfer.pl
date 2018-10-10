#!/usr/bin/perl
use strict;
use warnings;
open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<IN>){


	chomp;
	
	my @array=split(/\t/);
	for my $i(1..$#array) {
		print OUT "$array[0]\t$array[$i]\n";
		}
	
	}
close IN;
close OUT;
