# Requires SAMPLE_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

# list of files
files=( \
    "Aligned.toGenome.out.bam" \
    "Aligned.toTranscriptome.out.bam" \
    "feature_counts.txt" \
    "feature_counts.info" \
    "genes.fpkm_tracking" \
    "isoforms.fpkm_tracking" \
    "picard.CollectInsertSizes.txt" \
    "picard.MarkDuplicates.txt" \
    "picard.RnaSeqMetrics.txt" \
    "Signal.UniqueMultiple.str-.bw" \
    "Signal.UniqueMultiple.str+.bw" \
    "Signal.UniqueMultiple.both.bw" \
    "adapter_counts.info" \
    "rRNA_counts.info" \
    "trims.R1.fastq.gz" \
    "trims.R2.fastq.gz" \
    "kallisto_output/abundance.tsv" \
)

# list of sequins files
sequinsfiles=( \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_genes.tsv" \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv" \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info" \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats" \
    "anaquin_star/RnaAlign_summary.stats.info" \
)

# check files
for FILE in "${files[@]}"; do
    if [ ! -s $FILE ]; then
	echo "Missing $FILE"
	EXIT=1
    fi
done

# check sequins files
if [[ -n "$SEQUINS_REF" ]]; then
    for FILE in "${sequinsfiles[@]}"; do
        if [ ! -s $FILE ]; then
            echo "Missing $FILE"
            EXIT=1
        fi
    done
fi

if [[ $EXIT -ne 1 ]]; then
    python3 "$STAMPIPES/scripts/lims/upload_data.py" --aggregation_id ${AGGREGATION_ID} --complete_aggregation
fi

exit $EXIT
