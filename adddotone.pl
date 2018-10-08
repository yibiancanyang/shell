#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<>)
{
	$_=~ s/(>[A-Z0-9_]*\.1) .*?\n/$1/g;
	print OUT;
}
close IN;
close OUT;

