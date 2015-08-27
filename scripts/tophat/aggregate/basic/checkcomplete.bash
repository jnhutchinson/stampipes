# Requires LIBRARY_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

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
if [ ! -s $FILE ]; then
    echo "Missing $FILE"
    EXIT=1
fi
done

exit $EXIT
