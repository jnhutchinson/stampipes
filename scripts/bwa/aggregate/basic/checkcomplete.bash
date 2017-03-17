# Requires LIBRARY_NAME and GENOME to be in the environment 
# Checks that important files exist and are not size 0
# Checking for paired-end-specific files will only happen if PAIRED is set

EXIT=0

checkfile(){
    if [[ ! -s "$1" ]] ; then
        echo "Missing $1"
        EXIT=1
    fi
}

files=(
    "$LIBRARY_NAME.$GENOME.sorted.bam"
    "$LIBRARY_NAME.$GENOME.sorted.bam.bai"
    "$LIBRARY_NAME.tagcounts.txt" 
    "$LIBRARY_NAME.75_20.$GENOME.bw" 
    "$LIBRARY_NAME.75_20.$GENOME.uniques-density.bed.starch" 
    "$LIBRARY_NAME.75_20.normalized.$GENOME.bw" 
    "$LIBRARY_NAME.75_20.normalized.$GENOME.uniques-density.bed.starch" 
    "$LIBRARY_NAME.$GENOME.cuts.sorted.bed.starch" 
    "$LIBRARY_NAME.$GENOME.cutcounts.sorted.bed.starch" 
    "$LIBRARY_NAME.$GENOME.cutcounts.$READ_LENGTH.bw" 
    "$LIBRARY_NAME.MarkDuplicates.picard"
    )
paired_files=(
    "$LIBRARY_NAME.CollectInsertSizeMetrics.picard"
    "$LIBRARY_NAME.$GENOME.fragments.sorted.bed.starch"
    )

for FILE in "${files[@]}"; do
    checkfile "$FILE"
done

# Only check for InsertSizeMetrics on paired-end data
if [[ -n "$PAIRED" ]]; then
    for FILE in "${paired_files[@]}" ; do
        checkfile "$FILE"
    done
fi

if [[ $EXIT -ne 1 ]]; then
    python3 "$STAMPIPES/scripts/lims/upload_data.py" --aggregation_id ${AGGREGATION_ID} --complete_aggregation
fi

exit $EXIT
