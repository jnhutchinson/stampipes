# Checks that important files exist and are not size 0

EXIT=0

# list of files
files=( \
    "merged.transcriptome.cram" \
    "feature_counts.txt" \
    "genes.fpkm_tracking" \
    "isoforms.fpkm_tracking" \
    "kallisto.log" \
    "kallisto_adv.log" \
    "picard.MarkDuplicates.txt" \
    "picard.RnaSeqMetrics.txt" \
    "Signal.Unique.both.bw" \
    "Signal.Unique.str-.bw" \
    "Signal.Unique.str+.bw" \
    "adapter_counts.info" \
    "ribosomal_counts.info" \
    "kallisto_output/abundance.tsv" \
    "kallisto_output_adv/abundance.tsv" \
)

# Paired files only exist for paired-end aggregations.
paired_files=( \
    "picard.CollectInsertSizes.txt" \
)

# list of sequins files
# turned off until we get a sequins flag
sequins_files=( \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_genes.tsv" \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv" \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info" \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats" \
#    "anaquin_star/RnaAlign_summary.stats.info" \
)

function check_files() {
    for FILE in "$@" ; do
        if [ ! -s "$FILE" ]; then
            echo "Missing $FILE"
            EXIT=1
        fi
    done
}

# check files
check_files "${files[@]}"

if [[ -n "$PAIRED" ]] ; then
    check_files "${paired_files[@]}"
fi

if [[ -n "$SEQUINS_REF" ]]; then
    check_files "${sequins_files[@]}"
fi

if [[ $EXIT -eq 0 ]]; then
    python3 "$STAMPIPES/scripts/lims/upload_data.py" --aggregation_id "$AGGREGATION_ID" --complete_aggregation
fi

exit $EXIT
