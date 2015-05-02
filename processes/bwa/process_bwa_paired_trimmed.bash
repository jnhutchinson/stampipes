# Dependencies
source $MODULELOAD
module load bedops/2.4.2
module load bedtools/2.16.2
module load bwa/0.7.0
module load java/jdk1.7.0_05
module load picard/1.118
module load anaconda/2.1.0-2.7
module load samtools/1.2
module load gcc/4.7.2
module load R/3.1.0
module load git/2.3.3
module load coreutils/8.9

MAX_MISMATCHES=2
MIN_MAPPING_QUALITY=10
MAX_INSERT_SIZE=750

export FINAL_BAM=${SAMPLE_NAME}.sorted.bam
export UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam

export ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
export VERSION_FILE=${SAMPLE_NAME}.versions.txt

if [ -n "$REDO_ALIGNMENT" ]; then
    bash $STAMPIPES/scripts/bwa/reset_alignment.bash
fi

bash $STAMPIPES/scripts/versions.bash &> $VERSION_FILE

if [[ ( -n "$ADAPTER_P7" ) && ( -n "ADAPTER_P5" ) ]] ; then
  echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" > $ADAPTER_FILE
fi

if [ ! -e $FASTQ_DIR/${SAMPLE_NAME}_R1_fastqc -o ! -e $FASTQ_DIR/${SAMPLE_NAME}_R2_fastqc -o -n "$FORCE_FASTQ" ]; then
qsub -N ".fq${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  cd $FASTQ_DIR
  make -f $STAMPIPES/makefiles/fastqc.mk

  if [ "$UMI" = "True" ]; then
      echo "Tallying up top UMI tags seen in R1"
      zcat ${SAMPLE_NAME}_R1_???.fastq.gz | grep "^@" | cut -f 2 -d "+" | sort | uniq -c | sort -n -r | head -n 100 > ${SAMPLE_NAME}.topumis.txt
  fi

  bash /home/audrakj/stampipes/scripts/fastq/attachfiles.bash

  echo "FINISH: "
  date
__SCRIPT__
fi

NUMBER_FASTQ_FILES=`find $FASTQ_DIR -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l`
FASTQ_PAIR_HOLDS=""
FASTQ_PAIR_BAMS=""

# Indicate we have started this alignment and upload pertinent information
if [ ! -e ${FINAL_BAM} ]; then

python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
  -t ${LIMS_API_TOKEN} \
  --alignment_id ${ALIGNMENT_ID} \
  --start_alignment_progress \
  --adapter_file $ADAPTER_FILE \
  --version_file $VERSION_FILE

fi

for filenum in $(seq -f "%03g" 1 $NUMBER_FASTQ_FILES)
do
  NAME=".aln${SAMPLE_NAME}_${filenum}_${FLOWCELL}"
  BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"
    
if [ ! -e $BAMFILE -a ! -e ${FINAL_BAM} ]; then

qsub -l h_data=5650M -pe threads 1-2 -N ${NAME} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  if [ "$UMI" = "True" ]; then

  make -f $STAMPIPES/makefiles/bwa/bwa_paired_trimmed_umi.mk \
    FASTQ1_FILE=${FASTQ_DIR}/${SAMPLE_NAME}_R1_${filenum}.fastq.gz \
    FASTQ2_FILE=${FASTQ_DIR}/${SAMPLE_NAME}_R2_${filenum}.fastq.gz \
    OUTBAM=${BAMFILE} \
    TRIMSTATS=${SAMPLE_NAME}_${filenum}.trimstats.txt

  else

  make -f $STAMPIPES/makefiles/bwa/bwa_paired_trimmed.mk \
    FASTQ1_FILE=${FASTQ_DIR}/${SAMPLE_NAME}_R1_${filenum}.fastq.gz \
    FASTQ2_FILE=${FASTQ_DIR}/${SAMPLE_NAME}_R2_${filenum}.fastq.gz \
    OUTBAM=${BAMFILE} \
    TRIMSTATS=${SAMPLE_NAME}_${filenum}.trimstats.txt

  fi

  echo "FINISH: " `date`

__SCRIPT__

  # Only hold on alignments that are being run
  FASTQ_PAIR_HOLDS="$FASTQ_PAIR_HOLDS,$NAME"

fi

# Need to keep track of these even if they have already finshed
# for proper merging
FASTQ_PAIR_BAMS="${BAMFILE} ${FASTQ_PAIR_BAMS}"

done


if [ -n "$FASTQ_PAIR_HOLDS" ]; then
    SPLIT_ALIGN_HOLD="-hold_jid $FASTQ_PAIR_HOLDS"
fi

if [ ! -e ${FINAL_BAM} -o ! -e ${UNIQUES_BAM} ]; then

# If we are redoing this part, then we should make sure
# to redo all the other results as well
export FORCE_COUNTS=1
export REPROCESS=1

PROCESS_HOLD="-hold_jid .pb${SAMPLE_NAME}_${FLOWCELL}"

JOBNAME=".pb${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME"

# If we are processing UMI, we will need a lot of power for sorting!
if [ "$UMI" = "True" ]; then
  export SUBMIT_SLOTS="-pe threads 2"
fi

export PRE_DUP_BAM=${SAMPLE_NAME}.predup.sorted.bam
export PRE_DUP_BAM_PREFIX=${PRE_DUP_BAM%.*}
export FINAL_BAM_PREFIX=${FINAL_BAM%.*}

qsub ${SPLIT_ALIGN_HOLD} ${SUBMIT_SLOTS} -N "$JOBNAME" -V -cwd -S /bin/bash << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  if [ ! -e ${FINAL_BAM} ]; then
    echo "START MERGE: "
    date
  
    if [ "$NUMBER_FASTQ_FILES" -eq "1" ]; then
      mv ${SAMPLE_NAME}_001.sorted.bam ${FINAL_BAM}
    else
      samtools merge ${FINAL_BAM} ${FASTQ_PAIR_BAMS}
    fi
  
    # If we are working with a UMI, we can create the duplicate score and replace
    # the final BAM with one marked with duplicates
    if [ "$UMI" = "True" ]; then
      echo "START UMI DUP: "
      date
  
      rsync $FINAL_BAM \$TMPDIR/${PRE_DUP_BAM}
      rm $FINAL_BAM
  
      make -f $STAMPIPES/makefiles/umi/mark_duplicates.mk INPUT_BAM_FILE=\$TMPDIR/${PRE_DUP_BAM} \
        OUTPUT_BAM_FILE=\$TMPDIR/${FINAL_BAM}.presort
      samtools sort \$TMPDIR/${FINAL_BAM}.presort $FINAL_BAM_PREFIX
  
      echo "FINISH UMI DUP: "
      date
    fi

    echo "FINISH MERGE: "
    date
  else
    echo "$FINAL_BAM exists already"
  fi
  
  echo "START PROCESS BAM: "
  date

  make -f $STAMPIPES/makefiles/bwa/process_paired_bam.mk

  echo "FINISH PROCESS BAM: "
  date

  echo "START PICARD DUP: "
  date

  make -f $STAMPIPES/makefiles/picard/dups.mk

  echo "FINISH PICARD DUP: " date

  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --flowcell_lane_id ${FLOWCELL_LANE_ID} \
    --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard \
    --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard
  
  if [ "$NUMBER_FASTQ_FILES" -gt "1" ]
  then
    rm $FASTQ_PAIR_BAMS
  fi

  echo "FINISH: "
  date

__SCRIPT__

fi


if [ ! -e ${SAMPLE_NAME}.tagcounts.txt -o -n "$FORCE_COUNTS" ]; then

JOBNAME=".ct${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME"
   
qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: "
  hostname

  echo "START: " 
  date

  time bash $STAMPIPES/scripts/bwa/tagcounts.bash $SAMPLE_NAME $SAMPLE_NAME.sorted.bam $SAMPLE_NAME.tagcounts.txt $FASTQ_DIR
  # upload all data to the LIMS
  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --countsfile ${SAMPLE_NAME}.tagcounts.txt

  echo "FINISH: "
  date
      
__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out -o ! -e ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt ]; then

JOBNAME=".sp${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME" 
    
qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  make -f $STAMPIPES/makefiles/SPOT/spot-R1-paired.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
  # upload all data to the LIMS
  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt

  echo "FINISH: "
  date

__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.75_20.${GENOME}.bw ]; then

JOBNAME=".den${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME" 

qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  make -f $STAMPIPES/makefiles/densities/density.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME

  echo "FINISH: "
  date

__SCRIPT__

fi

if [ -n ${PROCESSING} ]; then

qsub -hold_jid ${PROCESSING} -N ".com${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: "
  hostname
  
  echo "START: "
  date

  bash $STAMPIPES/scripts/bwa/checkcomplete.bash

  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --finish_alignment

  bash /home/audrakj/stampipes/scripts/bwa/attachfiles.bash

  echo "FINISH: "
  date
__SCRIPT__

fi
