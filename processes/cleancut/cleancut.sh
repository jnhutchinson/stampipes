#!/bin/bash

set -e -u -o pipefail

module purge
module load jdk nextflow
module load anaconda
module load python/3.6.4
module load bcl2fastq2
module load bzip2
module load libpng

fcdir=$PWD

config=$STAMPIPES_DATA/cleancut/config.json
src="source.tsv"

if ! [[ -f "$config" ]] ; then
  echo "Config file '$config' does not exist"
  exit 1
fi

if ! [[ -f "$src" ]] ; then
  echo "Source file '$src' does not exist"
  exit 1
fi

python "$STAMPIPES/scripts/cleancut/create_cleancut_files.py"

# Use custom nextflow for conda support
nextflow-18.10 run "$STAMPIPES/scripts/cleancut/TalenTools/cleancut.nf" \
  -profile cluster,conda \
  -config "$STAMPIPES/processes/cleancut/nextflow.config" \
  --rundir "$fcdir" \
  --outdir "nextflow_output" \
  --manifest manifest.yaml \
  --specsheet spec-sheet.tsv \
  --minreads 0 \
  --config "$config" \
  -resume \
  "$@"
