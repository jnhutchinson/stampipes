# Requires LIBRARY_NAME and GENOME to be in the environment 
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "${LIBRARY_NAME}.${GENOME}.sorted.bam" \ 
    "${LIBRARY_NAME}.${GENOME}.sorted.bam.bai" \ 
    "${LIBRARY_NAME}.tagcounts.txt" \
    "${LIBRARY_NAME}.CollectInsertSizeMetrics.picard" \
    "${LIBRARY_NAME}.75_20.${GENOME}.bw" \
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
    python3 /home/audrakj/stampipes/scripts/lims/upload_data.py --aggregation_id 508 --complete_aggregation
fi

exit $EXIT
