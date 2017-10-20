# Uploads important files as attachments to objects on the LIMS
# Required environment names:
#  * STAMPIPES
#  * SAMPLE_NAME
#  * GENOME
#  * ALIGN_DIR
#  * ALIGNMENT_ID

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

ATTACH_ALIGNMENT="python3 $UPLOAD_SCRIPT --attach_file_contenttype SequencingData.flowcelllanealignment --attach_file_objectid ${ALIGNMENT_ID}"

cd ${ALIGN_DIR}

$ATTACH_ALIGNMENT --attach_directory ${ALIGN_DIR} --attach_file_purpose alignment-directory
$ATTACH_ALIGNMENT --attach_file Aligned.toTranscriptome.out.bam --attach_file_purpose all-alignments-bam  --attach_file_type bam
$ATTACH_ALIGNMENT --attach_file anaquin_star/RnaAlign_sequins.tsv --attach_file_purpose anaquin-rnaalign-star-sequins --attach_file_type plaintext
$ATTACH_ALIGNMENT --attach_file anaquin_star/RnaAlign_summary.stats --attach_file_purpose anaquin-rnaalign-star-summary --attach_file_type plaintext
$ATTACH_ALIGNMENT --attach_file trimmed/${SAMPLE_NAME}.trimmed.R1.fastq.gz --attach_file_purpose r1-fastq-trimmed --attach_file_type gzipped-fastq
$ATTACH_ALIGNMENT --attach_file trimmed/${SAMPLE_NAME}.trimmed.R2.fastq.gz --attach_file_purpose r2-fastq-trimmed --attach_file_type gzipped-fastq

#TODO: Attach more file types
