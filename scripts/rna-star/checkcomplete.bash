# Requires SAMPLE_NAME, GENOME and OUT_DIR to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "$OUT_DIR/${SAMPLE_NAME}.versions.txt" \
    "$OUT_DIR/Aligned.sortedByCoord.out.cram" \
    "$OUT_DIR/Aligned.toTranscriptome.out.cram" \
)

for FILE in "${files[@]}"; do
if [ ! -s $FILE ]; then
    echo "Missing $FILE"
    EXIT=1
fi
done

exit $EXIT
