# Requires LIBRARY_NAME, AGGREGATION_ID and GENOME to be in the environment
# Removes any of the listed files that exist

echo "RESETTING AGGREGATION ${AGGREGATION_ID} FOR ${LIBRARY_NAME}"

files=( \
    "${LIBRARY_NAME}.adapters.txt" \
    "${LIBRARY_NAME}.all.${GENOME}.bam" \
    "${LIBRARY_NAME}.all.${GENOME}.bam.bai" \
    "${LIBRARY_NAME}.all.${GENOME}.bw" \
    "${LIBRARY_NAME}.all.${GENOME}.starch" \
    "${LIBRARY_NAME}.neg.${GENOME}.bam" \
    "${LIBRARY_NAME}.neg.${GENOME}.bam.bai" \
    "${LIBRARY_NAME}.neg.${GENOME}.bw" \
    "${LIBRARY_NAME}.neg.${GENOME}.starch" \
    "${LIBRARY_NAME}.pos.${GENOME}.bam" \
    "${LIBRARY_NAME}.pos.${GENOME}.bam.bai" \
    "${LIBRARY_NAME}.pos.${GENOME}.bw" \
    "${LIBRARY_NAME}.pos.${GENOME}.starch" \
    "${LIBRARY_NAME}.bam_counts.txt" \
    "${LIBRARY_NAME}.read_counts.txt" \
    "${LIBRARY_NAME}.rRNAcounts.txt" \
    "${LIBRARY_NAME}.sample_summary.txt" \
    "${LIBRARY_NAME}.spotdups.txt" \
    "picard.${LIBRARY_NAME}.AlignmentSummary.txt" \
    "picard.${LIBRARY_NAME}.InsertSize.txt" \
    "picard.${LIBRARY_NAME}.RnaSeq.txt" \
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
