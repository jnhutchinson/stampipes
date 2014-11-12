#!/bin/bash
#TODO: Make this work
samplesheet=SampleSheet.csv
processing=processing.json

# jq scripts for reading $processing
JQ_GET_LANES=$(cat <<'__LANES__'
[.libraries[].lane]
| sort
| unique
| map(tostring)
| join (" ")
__LANES__
)
JQ_GET_BARCODES=$(cat <<'__JSON__'
.libraries[]
| select(.lane == ($lane | tonumber) )
| [.barcode_index, .samplesheet_name, .purpose ]
| join("\t")
__JSON__
)

# perl summation code
PERL_SUM=$(cat <<'__SUM__'
$f{$F[1]} += $F[0] ;
END {
    foreach $k (sort {$f{$b} <=> $f{$a} } keys %f ) {
        print "$k $f{$k}\n"
    }
}
__SUM__
)

# Okay let's actually process stuff now
lanes=$(jq  -r "$JQ_GET_LANES" "$processing" )

for lane in $lanes; do
    echo "============= Lane $lane ============="
    lanefiles=$(find Proj* Undet* -name "*L00${lane}_R1*.fastq.gz" | sed 's/.fastq.gz$/.barcodes.txt/')
    echo "Barcode,Clusters,Sample,Purpose"
    perl -ane "$PERL_SUM" $lanefiles \
    | sort \
    | join -a 1 - <( jq -r --arg lane "$lane" "$JQ_GET_BARCODES" "$processing" | sort ) \
    | sort -nrk2 \
    | awk '$3 || FNR < 30 {print} FNR==30 {print; print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"}' \
    | sed 's/\s\+/,/g; s/,/ /4g' # Restores spaces in "Purpose" field

    echo ''
done
