#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
while(<IN>)
{
	my $line = $_;
	$line =~ s/.*gene:(CGI_[0-9]*) .*\n/>$1\n/g;
	print OUT $line;
}
close IN;
close OUT;

