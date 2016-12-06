#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<>)
{
	$_=~ s/>.*\|ref\|(.*)\|.*\n/>$1\n/g;
	print OUT;
}
close IN;
close OUT;

