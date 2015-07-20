source $MODULELOAD
module load bedops/2.4.14
module load bowtie/1.0.0
module load cufflinks/2.2.1
module load gcc/4.7.2
module load java/jdk1.7.0_05
module load picard/1.118
module load samtools/1.2
module load tophat/2.0.13

source $PYTHON3_ACTIVATE
module load python/2.7.9

export SCRIPT_DIR="$STAMPIPES/scripts/tophat"
export REF_DIR=$(dirname "$BWAINDEX")
export ADAPTER_FILE="$SAMPLE_NAME.adapters.txt"
export VERSION_FILE="$SAMPLE_NAME.versions.txt"

filesize=$( du --total "$FASTQ_DIR"/*fastq.gz | tail -n1 | cut -f1)

# For small files, prioritize overall throughput
# For big files, we want them to finish someday.
if [ "$filesize" -gt 2000000 ] ; then
    THREADS=4
    SLOTS=2
    JOBS=4
else
    THREADS=1
    SLOTS=1
    JOBS=2
fi

bash $STAMPIPES/scripts/versions.bash &> $VERSION_FILE

qsub -cwd -V -q all.q -N .th${SAMPLE_NAME}_${FLOWCELL}_ALIGN#${ALIGNMENT_ID} -now no -pe threads $SLOTS -S /bin/bash <<__MAKE__
  set -x -e -o pipefail
  echo "Hostname: "
  hostname

  echo "START: "
  date

  if [[ ( -n "$ADAPTER_P7" ) && ( -n "$ADAPTER_P5" ) ]] ; then
    echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" > "$ADAPTER_FILE"
  fi

  python3 $STAMPIPES/scripts/lims/upload_data.py -a "$LIMS_API_URL" \
    -t "$LIMS_API_TOKEN" \
    --alignment_id "$ALIGNMENT_ID" \
    --start_alignment_progress \
    --adapter_file "$ADAPTER_FILE" \
    --version_file "$VERSION_FILE"

  make --keep-going -f \$STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk -j "$JOBS"

  bash $STAMPIPES/scripts/tophat/checkcomplete.bash

  bash $STAMPIPES/scripts/tophat/attachfiles.bash

  python3 $STAMPIPES/scripts/lims/upload_data.py -a "$LIMS_API_URL" \
    -t "$LIMS_API_TOKEN" \
    -f "$FLOWCELL" \
    --alignment_id "$ALIGNMENT_ID" \
    --finish_alignment

  echo "FINISH: "
  date

__MAKE__
