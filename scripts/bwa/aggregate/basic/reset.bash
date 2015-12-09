# Requires LIBRARY_NAME, AGGREGATION_ID and GENOME to be in the environment
# Removes any of the listed files that exist

echo "RESETTING AGGREGATION ${AGGREGATION_ID} FOR ${LIBRARY_NAME}"

files=( \
    "${LIBRARY_NAME}.${GENOME}.sorted.bam" \ 
    "${LIBRARY_NAME}.${GENOME}.sorted.bam.bai" \ 
    "${LIBRARY_NAME}.tagcounts.txt" \
    "${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam" \
    "${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam.bai" \
    "${LIBRARY_NAME}.CollectInsertSizeMetrics.picard" \
    "${LIBRARY_NAME}.CollectInsertSizeMetrics.picard.pdf" \
    "${LIBRARY_NAME}.MarkDuplicates.picard" \
    "${LIBRARY_NAME}.75_20.uniques-density.${GENOME}.bed.starch" \
    "${LIBRARY_NAME}.75_20.${GENOME}.bw" \
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.bw" \
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.uniques-density.bed.starch" \ 
    "$LIBRARY_NAME.${GENOME}.cuts.sorted.bed.starch" \
    "$LIBRARY_NAME.${GENOME}.cutcounts.sorted.bed.starch" \
    "$LIBRARY_NAME.${GENOME}.fragments.sorted.bed.starch" \
    "$LIBRARY_NAME.${GENOME}.cutcounts.$READ_LENGTH.bw" \
)

for FILE in "${files[@]}"; do
    if [ -e "$FILE" ]; then
        echo "Removing $FILE"
        rm $FILE
    fi
done

# Delete old job file logs

prefixes=( \
     ".AGG#${AGGREGATION_ID}" \
     ".agg" \
)

for PREFIX in "${prefixes[@]}"; do
    for FILE in `find . -name "$PREFIX${LIBRARY_NAME}*${FLOWCELL}.*"`; do
        echo "Removing $FILE"
        rm $FILE
    done
done

python3 $STAMPIPES/scripts/lims/upload_data.py --clear_aggregation_stats --aggregation_id ${AGGREGATION_ID}
