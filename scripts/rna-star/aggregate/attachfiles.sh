# Uploads important files as attachments to an aggregation object on the LIMS
# Required environment names:
#  * STAMPIPES
#  * LIBRARY_NAME
#  * GENOME
#  * AGGREGATION_DIR
#  * AGGREGATION_ID

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

ATTACH_AGGREGATION="python3 $UPLOAD_SCRIPT --attach_file_contenttype AggregationData.aggregation --attach_file_objectid ${AGGREGATION_ID}"

$ATTACH_AGGREGATION --attach_directory "$PWD/.." --attach_file_purpose aggregation-directory

# alignments
$ATTACH_AGGREGATION --attach_file merged.genome.bam --attach_file_purpose all-alignments-bam --attach_file_type bam
$ATTACH_AGGREGATION --attach_file merged.transcriptome.bam --attach_file_purpose transcriptome-alignments --attach_file_type bam

# cufflinks
$ATTACH_AGGREGATION --attach_file genes.fpkm_tracking --attach_file_purpose gene-expression-levels --attach_file_type plaintext
$ATTACH_AGGREGATION --attach_file isoforms.fpkm_tracking --attach_file_purpose isoform-expression-levels --attach_file_type plaintext

# kallisto
$ATTACH_AGGREGATION --attach_file kallisto_output/abundance.tsv --attach_file_purpose kallisto-abundance --attach_file_type tabtext

# featureCounts
$ATTACH_AGGREGATION --attach_file feature_counts.txt --attach_file_purpose featureCounts-output --attach_file_type plaintext

# densities
$ATTACH_AGGREGATION --attach_file Signal.UniqueMultiple.str+.bw --attach_file_purpose pos-coverage-bigwig --attach_file_type bigwig
$ATTACH_AGGREGATION --attach_file Signal.UniqueMultiple.str-.bw --attach_file_purpose neg-coverage-bigwig --attach_file_type bigwig
$ATTACH_AGGREGATION --attach_file Signal.UniqueMultiple.both.bw --attach_file_purpose all-coverage-bigwig --attach_file_type bigwig

# trimmed fastqs
#$ATTACH_AGGREGATION --attach_file trims.R1.fastq.gz --attach_file_purpose r1-fastq-trimmed --attach_file_type gzipped-fastq
#$ATTACH_AGGREGATION --attach_file trims.R2.fastq.gz --attach_file_purpose r2-fastq-trimmed --attach_file_type gzipped-fastq

# picard uploads
python3 $UPLOAD_SCRIPT --aggregation_id ${AGGREGATION_ID} --insertsfile picard.CollectInsertSizes.txt
python3 $UPLOAD_SCRIPT --aggregation_id ${AGGREGATION_ID} --dupsfile picard.MarkDuplicates.txt
