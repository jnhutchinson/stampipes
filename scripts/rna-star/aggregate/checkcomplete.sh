# Requires SAMPLE_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "Aligned.toGenome.out.bam" \
    "Aligned.toTranscriptome.out.bam" \
    "feature_counts.txt" \
    "genes.fpkm_tracking" \
    "isoforms.fpkm_tracking" \
    "picard.CollectInsertSizes.txt" \
    "picard.MarkDuplicates.txt" \
    "picard.RnaSeqMetrics.txt" \
    "Signal.UniqueMultiple.str-.bw" \
    "Signal.UniqueMultiple.str+.bw" \
    "Signal.UniqueMultiple.both.bw" \
    "adapter_counts.info" \
    "tagcounts.txt" \
    "trims.R1.fastq.gz" \
    "trims.R2.fastq.gz" \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_genes.tsv" \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv" \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info" \
    "anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats" \
    "anaquin_star/RnaAlign_summary.stats.info" \
    "kallisto_output/abundance.tsv" \
)

for FILE in "${files[@]}"; do
    if [ ! -s $FILE ]; then
	echo "Missing $FILE"
	EXIT=1
    fi
done

if [[ $EXIT -ne 1 ]]; then
    python3 "$STAMPIPES/scripts/lims/upload_data.py" --aggregation_id ${AGGREGATION_ID} --complete_aggregation
fi

exit $EXIT
