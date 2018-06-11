#!/bin/bash

set -e

here=$(dirname "$0")
cd "$here"

STAMPIPES=$(readlink -f ../..)
export STAMPIPES
nextflow run ../../processes/bwa/process_bwa_paired_trimmed.nf -profile test

# Verify

function cmp_picard() {
  name=$(basename "$1")
  expected=$(grep -v '^#' "expected/$name")
  actual=$(grep -v '^#' "output/$name")

  echo "Comparing $name..."
  diff <(echo "$expected") <(echo "$actual")
}

function cmp_starch() {
  name=$(basename "$1")
  if ! which unstarch &>/dev/null ; then
    echo "Cannot verify $name, unstarch is not available"
    return 0
  fi
  echo "Comparing $name..."
  cmp <(unstarch "expected/$name") <(unstarch "output/$name") \
    || (echo "$name does not match" ; false)
}

function cmp_bam() {
  name=$(basename "$1")
  if ! which samtools &>/dev/null ; then
    echo "Cannot verify $name, samtools is not available"
    return 0
  fi
  echo "Comparing $name..."
  cmp <(samtools view "expected/$name") <(samtools view "output/$name") \
    || (echo "$name does not match" ; false)

}

cmp_picard "MarkDuplicates.picard"
cmp_picard "CollectInsertSizeMetrics.picard"
cmp_picard "tagcounts.txt"
cmp_picard "subsample.spot.out"

cmp_starch "density.bed.starch"

cmp_bam "filtered.bam"
cmp_bam "marked.bam"
