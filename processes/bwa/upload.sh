#!/bin/bash

# Dependencies
source "$MODULELOAD"
module load zlib/1.2.8
module load bedops/2.4.35-typical
module load bedtools/2.25.0
module load bwa/0.7.12
module load jdk/1.8.0_92
module load picard/2.8.1
module load samtools/1.3
module load gcc/4.7.2
module load R/3.2.5
module load git/2.3.3
module load coreutils/8.25
module load pigz/2.3.3
module load modwt/1.0
module load htslib/1.6.0

# Load in this order specifically, currently the python3 activation
# overwrites the default "python" call, against advice
module load python/3.5.1
module load pysam/0.9.0
module load python/2.7.11

export MAX_MISMATCHES=2
export MIN_MAPPING_QUALITY=10
export MAX_INSERT_SIZE=750

JOB_BASENAME=${SAMPLE_NAME}_${FLOWCELL}_ALIGN#${ALIGNMENT_ID}

# expected output names
export FINAL_BAM=${SAMPLE_NAME}.sorted.bam
export FINAL_BAM_MARKED=${SAMPLE_NAME}.sorted.marked.bam # in case we ever want to maintain this...
export UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam
export ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
export VERSION_FILE=${SAMPLE_NAME}.versions.txt
export FASTQ_TMP=$ALIGN_DIR/fastq
export HIST=$SAMPLE_NAME.uniques.duphist.txt
export NUCLEAR_CHR=${NUCLEAR_CHR:-$BWAINDEX.nuclear.txt}

cd "$ALIGN_DIR"


python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --flowcell_lane_id ${FLOWCELL_LANE_ID} \
   --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard


python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --flowcell_lane_id ${FLOWCELL_LANE_ID} \
   --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard \
   --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard



# upload all data to the LIMS
python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --flowcell_lane_id ${FLOWCELL_LANE_ID} \
   --countsfile ${SAMPLE_NAME}.tagcounts.txt

   
if [[ -n "\$PAIRED" ]]; then
   python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt
else
   python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.rand.uniques.sorted.spotdups.txt
fi


	
python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --finish_alignment

bash $STAMPIPES/scripts/bwa/attachfiles.bash
