#!/bin/bash

inbam="$1"
outbam="$2"
numFrags="$3"
fragTotal="$4"

SUBSAMPLE_SEED=${SUBSAMPLE_SEED:-12345}

# A simpler version of random_sample.sh that discards read 2 and works quickly.

if [ -z "$fragTotal" ] ; then
    # Use samtools view to count the number of reads
    fragTotal=$(samtools view -c -F 128 "$inbam")
fi

if [[ "$numFrags" -ge "$fragTotal" ]] ; then
    samtools view -F 128 "$inbam" -o "$outbam"
    exit 0
fi

python "$STAMPIPES/scripts/bam/random_reads.py" \
    --singleend \
    --seed "$SUBSAMPLE_SEED" \
    <(samtools view -F 128 "$inbam" -u ) \
    "$outbam" \
    "$fragTotal" \
    "$numFrags"
