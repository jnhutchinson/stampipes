#!/bin/bash
module load jdk nextflow

if [[ $# -lt 2 ]] ; then
  echo "Usage: $0 readlengths /path/to/genome.fa [/path/to/genome2.fa ...]"
  echo
  echo "  readlengths: comma-separated list, like '36,76,101'"
  echo "  genomes:     full (absolute) paths to genome reference files"
  exit 1
fi

readlengths=$1
shift
genomes=$@

repo=$(readlink -f "$(dirname "$0")/../..")
GENOME_BUILD_DIR=${GENOME_BUILD_DIR:-$repo/genome_build}

cd "$GENOME_BUILD_DIR"

for g in "${genomes[@]}" ; do
  nuclear=${g/.fa/.nuclear.txt}
  nextflow run \
    "$repo/processes/build_references/build_reference.nf" \
    --genome "$g" \
    --outdir "$(basename "$g" .fa)" \
    --readlength "$readlengths" \
    -profile modules,cluster \
    -resume
done
