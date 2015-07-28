set -x -e -o pipefail

echo BAM FILES: $BAM_FILES

if [[ $BAM_COUNT -eq 1 ]]; then
    MERGETMP=${BAM_FILES}
else
    MERGETMP=${TMPDIR}/${LIBRARY_NAME}.mergetmp.sorted.bam
    samtools merge "${MERGETMP}" ${BAM_FILES}
fi

if [[ "$UMI" == "True" ]]; then
  echo "Using UMI mark dup"
  make -f $STAMPIPES/makefiles/umi/mark_duplicates.mk INPUT_BAM_FILE=${MERGETMP} OUTPUT_BAM_FILE=${FINAL_BAM} || rm ${FINAL_BAM}
else
  echo "Using Picard mark dup"
  make -f $STAMPIPES/makefiles/picard/dups.mk SAMPLE_NAME=${LIBRARY_NAME} BAMFILE="${MERGETMP}" OUTBAM=${FINAL_BAM} DUP_OUT=${LIBRARY_NAME}.MarkDuplicates.picard || rm ${FINAL_BAM}
fi
