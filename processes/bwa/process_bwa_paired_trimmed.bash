qsub -N "fq${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/fastqc.mk
__SCRIPT__

NUMBER_FASTQ_FILES=`find . -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l`
FASTQ_PAIR_HOLDS=""
FASTQ_PAIR_BAMS=""

for filenum in $(seq -f "%03g" 1 $NUMBER_FASTQ_FILES)
do
  NAME="aln${SAMPLE_NAME}_${filenum}"

qsub -l h_data=5650M -N ${NAME} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/bwa/bwa_paired_trimmed_make.mk \
    FASTQ1_FILE=${SAMPLE_NAME}_R1_${filenum}.fastq.gz \
    FASTQ2_FILE=${SAMPLE_NAME}_R2_${filenum}.fastq.gz \
    OUTBAM=${SAMPLE_NAME}_${filenum}.sorted.bam \
    TRIMSTATS=${SAMPLE_NAME}_${filenum}.trimstats.txt \
    ADAPTERFILE=$STAMPIPES/data/adapters/default.adapters
__SCRIPT__

  FASTQ_PAIR_HOLDS="$FASTQ_PAIR_HOLDS $NAME"
  FASTQ_PAIR_BAMS="${SAMPLE_NAME}_${filenum}.sorted.bam ${FASTQ_PAIR_BAMS}"
done

qsub -hold_jid $FASTQ_PAIR_HOLDS -N "pb${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  if [ "$NUMBER_FASTQ_FILES" -eq "1" ]
  then
    mv ${SAMPLE_NAME}_001.sorted.bam ${SAMPLE_NAME}.sorted.bam
  else
    samtools merge ${SAMPLE_NAME}.sorted.bam ${FASTQ_PAIR_BAMS}
  fi
  
  make -f $STAMPIPES/makefiles/bwa/process_paired_bam.mk
  
  if [ "$NUMBER_FASTQ_FILES" -gt "1" ]
  then
    rm $FASTQ_PAIR_BAMS
  fi
__SCRIPT__

qsub -hold_jid "pb${SAMPLE_NAME}" -N "ct${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  bash $STAMPIPES/scripts/bwa/tagcounts.bash $SAMPLE_NAME $SAMPLE_NAME.sorted.bam $SAMPLE_NAME.tagcounts.txt
  # upload all data to the LIMS
  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --countsfile ${SAMPLE_NAME}.tagcounts.txt \
      --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard
      
__SCRIPT__

qsub -hold_jid "pb${SAMPLE_NAME}" -N "sp${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  make -f $STAMPIPES/makefiles/SPOT/spot-R1.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
  #make -f $STAMPIPES/makefiles/SPOT/spot-Rboth.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
  # upload all data to the LIMS
  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
      --dupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt
__SCRIPT__

qsub -hold_jid "pb${SAMPLE_NAME}" -N "den${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  make -f $STAMPIPES/makefiles/densities/density.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
__SCRIPT__
