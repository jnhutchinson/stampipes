#!/usr/bin/perl

use strict;
use warnings;

my %starts = ();
my %ends = ();
my %changes = ();

my $curChrom = "";
my $farthestEnd = 0;


sub catchup {
  my @pos =  sort { $a <=> $b } grep { $changes{$_} != 0 } keys %changes ;
  my $density = 0;
  for (my $i = 0; $i < @pos - 1; $i++ ){
    $density += $changes{$pos[$i]};
    if ($density > 0){
      print "$curChrom\t$pos[$i]\t$pos[$i+1]\t$density\n";
    }
  }
}

while (<>){
  chomp;
  my ($chrom, $start, $end) = split;

  if ($start > $farthestEnd || $curChrom ne $chrom ) {
    catchup();
    %changes = ();
  }
  $changes{$start}++;
  $changes{$end}--;

  $farthestEnd = ($end > $farthestEnd) ?  $end : $farthestEnd;
  $curChrom = $chrom;

}

catchup();
