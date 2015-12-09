# Requires LIBRARY_NAME and GENOME to be in the environment 
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "${LIBRARY_NAME}.${GENOME}.sorted.bam" \ 
    "${LIBRARY_NAME}.${GENOME}.sorted.bam.bai" \ 
    "${LIBRARY_NAME}.tagcounts.txt" \
    "${LIBRARY_NAME}.CollectInsertSizeMetrics.picard" \
    "${LIBRARY_NAME}.75_20.${GENOME}.bw" \
    "${LIBRARY_NAME}.75_20.${GENOME}.uniques-density.bed.starch" \
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.bw" \
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.uniques-density.bed.starch" \
    "$LIBRARY_NAME.${GENOME}.cuts.sorted.bed.starch" \
    "$LIBRARY_NAME.${GENOME}.cutcounts.sorted.bed.starch" \
    "$LIBRARY_NAME.${GENOME}.fragments.sorted.bed.starch" \
    "$LIBRARY_NAME.${GENOME}.cutcounts.$READ_LENGTH.bw" \
)

for FILE in "${files[@]}"; do
if [ ! -s $FILE ]; then
    echo "Missing $FILE"
    EXIT=1
fi
done

# Only check for MarkDuplicates report if this is not a UMI sample
if [[ "$UMI" != "True" ]]; then
    if [ ! -s "${LIBRARY_NAME}.MarkDuplicates.picard" ]; then
        echo "Missing $FILE"
        EXIT=1
    fi
fi

if [[ $EXIT -ne 1 ]]; then
    python3 "$STAMPIPES/scripts/lims/upload_data.py" --aggregation_id ${AGGREGATION_ID} --complete_aggregation
fi

exit $EXIT
