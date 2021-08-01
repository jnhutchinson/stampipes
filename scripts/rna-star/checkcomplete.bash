# Requires SAMPLE_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "${SAMPLE_NAME}.versions.txt" \
    "output/Aligned.sortedByCoord.out.bam" \
    "output/Aligned.toTranscriptome.out.bam" \
)

for FILE in "${files[@]}"; do
if [ ! -s $FILE ]; then
    echo "Missing $FILE"
    EXIT=1
fi
done

exit $EXIT
