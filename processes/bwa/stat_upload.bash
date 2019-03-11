  cd $ALIGN_DIR
export FINAL_BAM=${SAMPLE_NAME}.sorted.bam
export UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam

export ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
export VERSION_FILE=${SAMPLE_NAME}.versions.txt
  bash $STAMPIPES/scripts/bwa/checkcomplete.bash
  
  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt  
  
  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --flowcell_lane_id ${FLOWCELL_LANE_ID} \
    --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard \
    --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard  
  
  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --finish_alignment

  bash $STAMPIPES/scripts/bwa/attachfiles.bash

    python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --countsfile ${SAMPLE_NAME}.tagcounts.txt
      
