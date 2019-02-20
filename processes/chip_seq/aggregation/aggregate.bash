#!/bin/bash

version=2.0.1

cd "$(dirname "$0")"

source "$MODULELOAD"
module purge
module load jdk
module load nextflow
module load python/3.5.1

if [[ $(wc -w <<< "$BAM_FILES") -gt 1 ]] ; then
  bamfiles="$(sed 's/\s\+/,/g' <<< "$BAM_FILES")"
else
  bamfiles=$BAM_FILES
fi
outdir="output_${version}"
workdir="work"

WIN=75
BINI=20


export MAPPABLE_REGIONS=${MAPPABLE_REGIONS:-$GENOME_INDEX.K${READ_LENGTH}.mappable_only.bed}
export NUCLEAR_CHR=${NUCLEAR_CHR:-$GENOME_INDEX.nuclear.txt}

# Remove old stuff if necessary
if [[ -n "$REDO_AGGREGATION" ]] ; then
  rm -rf "$outdir"
  rm -rf "$workdir"
  echo \
  python3 "$STAMPIPES/scripts/lims/upload_data.py" \
    --clear_aggregation_stats \
    --aggregation_id "$AGGREGATION_ID"
fi

# Tell LIMS we're starting alignment
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL" \
  -t "$LIMS_API_TOKEN" \
  --aggregation_id "$AGGREGATION_ID" \
  --start_aggregation

# Run the whole process
nextflow run \
  "$STAMPIPES/processes/chip_seq/aggregation/aggregate.nf" \
  -c "$STAMPIPES/nextflow.config" \
  -w "$workdir" \
  --treatment "$bamfiles" \
  --genome "$GENOME_INDEX" \
  --mappable "$MAPPABLE_REGIONS" \
  --UMI "$UMI" \
  --readlength "$READ_LENGTH" \
  --outdir "$outdir" \
  --threads 3 \
  -profile cluster,modules \
  -with-report nextflow.report.html \
  -with-dag nextflow.flowchart.html \
  -with-timeline nextflow.timeline.html \
  -resume

#( cd "$outdir" \
#    && bash "$STAMPIPES/scripts/chip_seq/aggregate/basic/attachfiles_nextflow.bash" \
#    && bash "$STAMPIPES/scripts/chip_seq/aggregate/basic/uploadcounts_nextflow.bash"
#)
# Upload results
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL" \
  -t "$LIMS_API_TOKEN" \
  --aggregation_id "$AGGREGATION_ID" \
  --complete_aggregation

