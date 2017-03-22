source $MODULELOAD
module load bedops/2.4.19
module load jdk/1.8.0_92
module load gcc/4.7.2
module load R/3.2.5
module load picard/2.8.1
module load samtools/1.3
module load git/2.3.3
module load coreutils/8.25
module load modwt/1.0
module load hotspot2/20161006

module load python/3.5.1
module load pysam/0.9.0

WIN=75
BINI=20

cd $AGGREGATION_FOLDER
BAM_COUNT=`ls $BAM_FILES | wc -l`
QUEUE=queue0

JOB_BASENAME=".AGG#${AGGREGATION_ID}"

export FINAL_BAM=${LIBRARY_NAME}.${GENOME}.sorted.bam
export FINAL_BAM_MARKED=${LIBRARY_NAME}.${GENOME}.sorted.marked.bam
export FINAL_UNIQUES_BAM=${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam
export TAGCOUNTS_FILE=${LIBRARY_NAME}.tagcounts.txt
export DENSITY_STARCH=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.uniques-density.bed.starch
export DENSITY_BIGWIG=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.bw
export NORM_DENSITY_BIGWIG=${LIBRARY_NAME}.${WIN}_${BINI}.normalized.${GENOME}.bw
export CUTCOUNTS_BIGWIG=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cutcounts.$READ_LENGTH.bw
export INSERT_FILE=${LIBRARY_NAME}.CollectInsertSizeMetrics.picard
export DUPS_FILE=${LIBRARY_NAME}.MarkDuplicates.picard

export HOTSPOT_DIR=peaks
HOTSPOT_PREFIX=$(basename "$FINAL_UNIQUES_BAM" .bam)
export HOTSPOT_CALLS=$HOTSPOT_DIR/$HOTSPOT_PREFIX.hotspots.fdr0.05.starch
export HOTSPOT_DENSITY=$HOTSPOT_DIR/$HOTSPOT_PREFIX.sorted.density.bw

export HOTSPOT_SCRIPT="hotspot2.sh"
export MAPPABLE_REGIONS=${MAPPABLE_REGIONS:-$GENOME_INDEX.K${READ_LENGTH}.mappable_only.bed}
export CHROM_SIZES=${CHROM_SIZES:-$GENOME_INDEX.chrom_sizes.bed}
export CENTER_SITES=${CENTER_SITES:-$GENOME_INDEX.K${READ_LENGTH}.center_sites.n100.starch}

if [ -n "$REDO_AGGREGATION" ]; then
    bash $STAMPIPES/scripts/bwa/aggregate/basic/reset.bash
    rm -r $HOTSPOT_DIR
fi

MERGE_DUP_JOBNAME=${JOB_BASENAME}_merge_dup
COUNT_JOBNAME=${JOB_BASENAME}_count
HOTSPOT_JOBNAME=${JOB_BASENAME}_hotspot
DENSITY_JOBNAME=${JOB_BASENAME}_density
CUTCOUNTS_JOBNAME=${JOB_BASENAME}_cutcounts

# Check out files match first
#python3 $STAMPIPES/scripts/utility/md5check.py bamfiles.txt || exit 1

PROCESSING=""

if [[ ! -s "$FINAL_BAM.bai" || ! -s "$FINAL_UNIQUES_BAM.bai" || ! -s "$DUPS_FILE" || (! -s "$INSERT_FILE" && -n "$PAIRED") ]]; then

PROCESSING="$PROCESSING,${MERGE_DUP_JOBNAME}"

qsub ${SUBMIT_SLOTS} -N "${MERGE_DUP_JOBNAME}" -q ${QUEUE} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: processing BAMs"
  date

  # merge BAMs
  if [[ ! -s "$FINAL_BAM" ]] ; then
    
    if [[ $BAM_COUNT -eq 1 ]]; then
      rsync ${BAM_FILES} ${FINAL_BAM}
    else
      samtools merge ${FINAL_BAM} ${BAM_FILES}
    fi
    samtools index ${FINAL_BAM}

  fi

  echo "PROCESSING: merged BAMs"
  date

  # process BAM
  if [[ ! -s $FINAL_UNIQUES_BAM ]] ; then

    if [[ "$UMI" == "True" ]]; then
      make -f "$STAMPIPES/makefiles/picard/dups_cigarumi.mk" SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
      mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
      samtools view -b -F 1536 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
    else
      make -f "$STAMPIPES/makefiles/picard/dups_cigar.mk" SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
      mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
      samtools view -b -F 512 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
    fi

    samtools index ${FINAL_BAM}
    samtools index ${FINAL_UNIQUES_BAM}

  fi

  echo "PROCESSING: marked duplicates"
  date
  
  # calculate insert sizes
  if [[ ! -s $INSERT_FILE ]]; then
    make -f "$STAMPIPES/makefiles/picard/insert_size_metrics.mk" "SAMPLE_NAME=${LIBRARY_NAME}" "BAMFILE=${FINAL_UNIQUES_BAM}" "INSERTMETRICS=${INSERT_FILE}"
  fi

  # placeholder for QC metric script

  echo "FINISH: processing BAMs"
  date

__SCRIPT__

fi

# Run Hotspot2
if [[ ! -s "$HOTSPOT_CALLS" || ! -s "$HOTSPOT_DENSITY" ]] ; then
  PROCESSING="$PROCESSING,${HOTSPOT_JOBNAME}"
  HOTSPOT_SPOT=$HOTSPOT_DIR/$HOTSPOT_PREFIX.SPOT.txt
  qsub ${SUBMIT_SLOTS} -hold_jid "${MERGE_DUP_JOBNAME}" -N "${HOTSPOT_JOBNAME}" -q ${QUEUE} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
    "$HOTSPOT_SCRIPT"  -F 0.5 -s 12345 -M "$MAPPABLE_REGIONS" -c "$CHROM_SIZES" -C "$CENTER_SITES" "$FINAL_UNIQUES_BAM"  "$HOTSPOT_DIR"
    "$STAMPIPES/scripts/SPOT/info.sh" "$HOTSPOT_CALLS" hotspot2 \$(cat $HOTSPOT_SPOT) > "$HOTSPOT_PREFIX.hotspot2.info"
__SCRIPT__
fi

if [ ! -e $TAGCOUNTS_FILE ]; then

PROCESSING="$PROCESSING,${COUNT_JOBNAME}"

qsub -N "${COUNT_JOBNAME}" -q ${QUEUE} -hold_jid ${MERGE_DUP_JOBNAME} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  python3 $STAMPIPES/scripts/bwa/bamcounts.py $FINAL_BAM $TAGCOUNTS_FILE

  echo "FINISH: "
  date

__SCRIPT__

fi

if [[ ! -s "$DENSITY_BIGWIG" || ! -s "$NORM_DENSITY_BIGWIG" ]]; then

PROCESSING="$PROCESSING,${DENSITY_JOBNAME}"

qsub -N "${DENSITY_JOBNAME}" -q ${QUEUE} -hold_jid "${MERGE_DUP_JOBNAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  make -f $STAMPIPES/makefiles/densities/density.mk FAI=${GENOME_INDEX}.fai SAMPLE_NAME=${LIBRARY_NAME} GENOME=${GENOME} \
     BAMFILE=${FINAL_UNIQUES_BAM} STARCH_OUT=${DENSITY_STARCH} BIGWIG_OUT=${DENSITY_BIGWIG}
  make -f "$STAMPIPES/makefiles/densities/normalize-density.mk" BAMFILE=${FINAL_UNIQUES_BAM} SAMPLE_NAME=${LIBRARY_NAME} FAI=${GENOME_INDEX}.fai

  echo "FINISH: "
  date

__SCRIPT__

fi

if [ ! -e $CUTCOUNTS_BIGWIG ]; then

PROCESSING="$PROCESSING,${CUTCOUNTS_JOBNAME}"

qsub -N "${CUTCOUNTS_JOBNAME}" -q ${QUEUE} -hold_jid "${MERGE_DUP_JOBNAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__

  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  bash $STAMPIPES/scripts/bwa/aggregate/basic/cutcounts.bash

  echo "FINISH: "
  date

__SCRIPT__

fi

if [ -n ${PROCESSING} ]; then

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py
ATTACH_AGGREGATION="python3 $UPLOAD_SCRIPT --attach_file_contenttype AggregationData.aggregation --attach_file_objectid ${AGGREGATION_ID}"
qsub -N "${JOB_BASENAME}_complete" -q ${QUEUE} -hold_jid "${PROCESSING}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  bash $STAMPIPES/scripts/bwa/aggregate/basic/checkcomplete.bash
  bash $STAMPIPES/scripts/bwa/aggregate/basic/attachfiles.bash

  echo "FINISH: "
  date

__SCRIPT__

fi
