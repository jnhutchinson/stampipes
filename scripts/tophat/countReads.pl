#!/usr/bin/perl -w
use strict;

# add up lines in fastq files, print name key value

if (@ARGV < 3) {
    die "Usage:  $0 name key <readfile> [readfiles...]\n";
}

my ($name,$key,@readFiles) = @ARGV;

my $sum = 0;
foreach my $file (@readFiles) {
    chomp( my $lines = `zcat $file | wc -l` );
    $sum += ($lines / 4);
}

print "$name\t$key\t$sum\n";

