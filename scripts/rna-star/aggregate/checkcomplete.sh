# Requires SAMPLE_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "Aligned.toGenome.out.bam" \
    "Aligned.toTranscriptome.out.bam" \
    "feature_counts.txt" \
    "genes.fpkm_tracking" \
    "isoforms.fpkm_tracking" \
    "Signal.UniqueMultiple.str-.bw" \
    "Signal.UniqueMultiple.str+.bw" \
    "Signal.UniqueMultiple.both.bw" \
    "trims.R1.fastq.gz" \
    "trims.R2.fastq.gz" \
    "anaquin_cufflinks/RnaExpression_genes.tsv" \
    "anaquin_cufflinks/RnaExpression_isoforms.tsv" \
    "anaquin_cufflinks/RnaExpression_summary.stats" \
    "anaquin_kallisto/RnaExpression_genes.tsv" \
    "anaquin_kallisto/RnaExpression_isoforms.tsv" \
    "anaquin_kallisto/RnaExpression_summary.stats" \
    "anaquin_star/RnaAlign_summary.stats" \
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
