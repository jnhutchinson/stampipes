#!/bin/bash
#TODO: Make this work
samplesheet=SampleSheet.csv
lanes="L001"
for lane in $lanes; do
    echo "Lane $lane"
    echo "========="
    lanefiles=$(find Proj* Undet* -name "*${lane}_R1*.fastq.gz" | sed 's/.fastq.gz$/.barcodes.txt/')
    echo -e "Barcode,Clusters,Sample_Assigned"
    perl -ane '$f{$F[0]} += $F[1] ; END { foreach $k (sort {$f{$b} <=> $f{$a} } keys %f ) { print "$k $f{$k}\n" }}' $lanefiles \
    | sort \
    | join -a 1 - <( awk -v LANE=$lane 'BEGIN{ FS=","; OFS="	"} ("L00" $2 == LANE ) {print $5, $3 }' $samplesheet  | sort ) \
    | sort -nrk2 \
    | awk '$3 || FNR < 30 {print} FNR==30 {print; print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"}' \
    | sed 's/\s\+/,/g'

    echo ''
done
