#!/usr/bin/env bash
# 
USAGE_MSG="Usage: $0 in.bam out.bam numPairs [pairTotal]"

source $MODULELOAD
module load samtools/1.2
source $PYTHON3_ACTIVATE

die(){
  echo $@ >&2
  exit 1
}

if [ $# -lt 3 ] ; then
  die "$USAGE"
fi

inbam="$1"
outbam="$2"
numPairs="$3"
pairTotal="$4"

# If not specified... calculate from the index file
if [ -z "$pairTotal" ] ; then
  if [ ! -s "$inbam.bai" ] ;then
    samtools index $inbam
  fi

  readTotal=$(samtools idxstats $inbam | awk '{sum+=$3}END{print sum}')
  pairTotal=$((readTotal / 2))
fi

echo "Random sampling started at: `date +%D%t%T.%N`"
echo "Sampling $numPairs read-pairs from $INBAMFILE to $OUTBAMFILE (from $pairTotal total)..."

python3 $STAMPIPES/scripts/bam/random-reads.py "$inbam" "$outbam" "$pairTotal" "$numPairs"

echo "Random sampling finished at: `date +%D%t%T.%N`"
