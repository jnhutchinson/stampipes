# Requires SAMPLE_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

# list of files
files=( \
    "Aligned.toGenome.out.bam" \
    "Aligned.toGenome.out.bam.bai" \
    "Aligned.toTranscriptome.out.bam" \
    "feature_counts.txt" \
    "feature_counts.txt.summary" \
    "genes.fpkm_tracking" \
    "isoforms.fpkm_tracking" \
    "kallisto.log" \
    "kallisto_adv.log" \
    "picard.CollectInsertSizes.txt" \
    "picard.MarkDuplicates.txt" \
    "picard.RnaSeqMetrics.txt" \
    "Signal.Unique.both.bw" \
    "Signal.Unique.str-.bw" \
    "Signal.Unique.str+.bw" \
    "adapter_counts.info" \
    "ribosomal_counts.info" \
    "trims.R1.fastq.gz" \
    "trims.R2.fastq.gz" \
    "kallisto_output/abundance.tsv" \
    "kallisto_output_adv/abundance.tsv" \
)

# list of sequins files
# turned off until we get a sequins flag
#sequinsfiles=( \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_genes.tsv" \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv" \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info" \
#    "anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats" \
#    "anaquin_star/RnaAlign_summary.stats.info" \
#)

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
