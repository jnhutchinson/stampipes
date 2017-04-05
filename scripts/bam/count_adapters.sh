#!/bin/bash

module load bwa
module load samtools
module load jdk
module load picard/2.9.0

bam=${1:-/dev/stdin}
adapters=${2:-$STAMPIPES/data/adapters/alladapters.fa}

if [ ! -e "$adapters.bwt" ] ; then
  bwa index "$adapters"
fi

samtools view -u -f4 "$bam" \
  | samtools fastq - \
  | picard FifoBuffer \
  | bwa mem -t 2 -k 13 "$adapters" - \
  | samtools view -F4 -c
