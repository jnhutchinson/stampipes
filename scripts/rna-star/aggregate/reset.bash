# Requires LIBRARY_NAME, AGGREGATION_ID and GENOME to be in the environment
# Removes any of the listed files that exist

echo "RESETTING AGGREGATION ${AGGREGATION_ID} FOR ${LIBRARY_NAME}"

files=( \
"Aligned.toGenome.out.bam" \
"Aligned.toGenome.out.bam.bai" \
"Aligned.toTranscriptome.out.bam" \
"feature_counts.txt" \
"feature_counts.txt.summary" \
"genes.fpkm_tracking" \
"isoforms.fpkm_tracking" \
"Signal.UniqueMultiple.str-.bw" \
"Signal.UniqueMultiple.str+.bw" \
"Signal.UniqueMultiple.both.bw" \
"Signal.UniqueMultiple.str-.starch" \
"Signal.UniqueMultiple.str+.starch" \
"Signal.UniqueMultiple.both.starch" \
"Signal.Unique.str-.bw" \
"Signal.Unique.str+.bw" \
"Signal.Unique.both.bw" \
"Signal.Unique.str-.starch" \
"Signal.Unique.str+.starch" \
"Signal.Unique.both.starch" \
"skipped.gtf" \
"transcripts.gtf" \
"trims.R1.fastq.gz" \
"trims.R2.fastq.gz" \
)

dirs=( \
    "anaquin_cufflinks" \
    "anaquin_kallisto"  \
    "anaquin_star"      \
    "kallisto_output"   \
)

for FILE in "${files[@]}"; do
    if [ -e "$FILE" ]; then
        echo "Removing $FILE"
        rm $FILE
    fi
done

for DIR in "${files[@]}"; do
    if [ -d $DIR ]; then
	echo "Removing Directory $DIR"
        rm -r $DIR
    fi
done

python3 $STAMPIPES/scripts/lims/upload_data.py --clear_aggregation_stats --aggregation_id ${AGGREGATION_ID}
