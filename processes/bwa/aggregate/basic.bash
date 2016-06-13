
source $MODULELOAD
module load bedops/2.4.15
module load java/jdk1.7.0_05
module load gcc/4.7.2
module load R/3.1.0
module load picard/1.118
module load samtools/1.2
module load git/2.3.3
module load coreutils/8.9

WIN=75
BINI=20

cd $AGGREGATION_FOLDER

BAM_COUNT=`ls $BAM_FILES | wc -l`

JOB_BASENAME=".AGG#${AGGREGATION_ID}"

export FINAL_BAM=${LIBRARY_NAME}.${GENOME}.sorted.bam
export TEMP_UNIQUES_BAM=${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam
export TAGCOUNTS_FILE=${LIBRARY_NAME}.tagcounts.txt
export DENSITY_STARCH=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.uniques-density.bed.starch
export DENSITY_BIGWIG=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.bw
export CUTCOUNTS_BIGWIG=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cutcounts.$READ_LENGTH.bw

if [ -n "$REDO_AGGREGATION" ]; then
    bash $STAMPIPES/scripts/bwa/aggregate/basic/reset.bash
fi

MERGE_DUP_JOBNAME=${JOB_BASENAME}_merge_dup
COUNT_JOBNAME=${JOB_BASENAME}_count
DENSITY_JOBNAME=${JOB_BASENAME}_density
CUTCOUNTS_JOBNAME=${JOB_BASENAME}_cutcounts

# Check out files match first
#python3 $STAMPIPES/scripts/utility/md5check.py bamfiles.txt || exit 1

PROCESSING=""

if [ ! -e ${FINAL_BAM} ]; then

PROCESSING="$PROCESSING,${MERGE_DUP_JOBNAME}"

# If we are  UMI, we will need a lot of power for sorting!
if [ "$UMI" = "True" ]; then
  echo "Processing with UMI"
  export SUBMIT_SLOTS="-pe threads 4"
  echo "Excluding duplicates from uniques bam"
  export EXCLUDE_FLAG=1536
fi

qsub ${SUBMIT_SLOTS} -N "${MERGE_DUP_JOBNAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  export BAM_COUNT=$BAM_COUNT

  bash $STAMPIPES/scripts/bwa/aggregate/merge.bash

  make -f $STAMPIPES/makefiles/bwa/process_paired_bam.mk SAMPLE_NAME=${LIBRARY_NAME} INBAM=${FINAL_BAM} OUTBAM=${TEMP_UNIQUES_BAM} info uniques

  echo "FINISH: "
  date

__SCRIPT__

fi

#if [ ! -e $TAGCOUNTS_FILE ]; then
#
#PROCESSING="$PROCESSING,${COUNT_JOBNAME}"
#
#qsub -N "${COUNT_JOBNAME}" -hold_jid ${MERGE_DUP_JOBNAME} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
#  set -x -e -o pipefail
#
#  echo "Hostname: "
#  hostname
#
#  echo "START: "
#  date
#
#  python3 $STAMPIPES/scripts/bwa/bamcounts.py $FINAL_BAM $TAGCOUNTS_FILE
#
#  echo "FINISH: "
#  date
#
#__SCRIPT__
#
#fi
#
#if [ ! -e $DENSITY_BIGWIG ]; then
#
#PROCESSING="$PROCESSING,${DENSITY_JOBNAME}"
#
#qsub -N "${DENSITY_JOBNAME}" -hold_jid "${MERGE_DUP_JOBNAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
#  set -x -e -o pipefail
#
#  echo "Hostname: "
#  hostname
#
#  echo "START: "
#  date
#
#  make -f $STAMPIPES/makefiles/densities/density.mk FAI=${GENOME_INDEX}.fai SAMPLE_NAME=${LIBRARY_NAME} GENOME=${GENOME} \
#     BAMFILE=${TEMP_UNIQUES_BAM} STARCH_OUT=${DENSITY_STARCH} BIGWIG_OUT=${DENSITY_BIGWIG}
#  make -f "$STAMPIPES/makefiles/densities/normalize-density.mk" BAMFILE=${TEMP_UNIQUES_BAM} SAMPLE_NAME=${LIBRARY_NAME} FAI=${GENOME_INDEX}.fai
#
#  echo "FINISH: "
#  date
#
#__SCRIPT__
#
#fi
#
#if [ ! -e $CUTCOUNTS_BIGWIG ]; then
#
#PROCESSING="$PROCESSING,${CUTCOUNTS_JOBNAME}"
#
#qsub -N "${CUTCOUNTS_JOBNAME}" -hold_jid "${MERGE_DUP_JOBNAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
#
#  set -x -e -o pipefail
#
#  echo "Hostname: "
#  hostname
#
#  echo "START: "
#  date
#
#  bash $STAMPIPES/scripts/bwa/aggregate/basic/cutcounts.bash 
#
#  echo "FINISH: "
#  date
#
#__SCRIPT__
#
#fi

if [ -n ${PROCESSING} ]; then

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py
ATTACH_AGGREGATION="python3 $UPLOAD_SCRIPT --attach_file_contenttype AggregationData.aggregation --attach_file_objectid ${AGGREGATION_ID}"
qsub -N "${JOB_BASENAME}_complete" -hold_jid "${PROCESSING}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  #bash $STAMPIPES/scripts/bwa/aggregate/basic/checkcomplete.bash
  #bash $STAMPIPES/scripts/bwa/aggregate/basic/attachfiles.bash
  $ATTACH_AGGREGATION --attach_directory `pwd` --attach_file_purpose aggregation-directory 
  $ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.${GENOME}.sorted.bam --attach_file_purpose all-alignments-bam --attach_file_type bam
  $ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.${GENOME}.sorted.bam.bai --attach_file_purpose bam-index --attach_file_type bam-index

  #rm ${TEMP_UNIQUES_BAM} ${TEMP_UNIQUES_BAM}.bai

  echo "FINISH: "
  date

__SCRIPT__

fi
