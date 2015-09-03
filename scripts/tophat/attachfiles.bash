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
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.all.$GENOME.bam      --attach_file_purpose all-alignments-bam  --attach_file_type bam
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.all.$GENOME.bam.bai  --attach_file_purpose bam-index           --attach_file_type bam-index
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.pos.$GENOME.bam      --attach_file_purpose pos-alignments-bam  --attach_file_type bam
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.pos.$GENOME.bam.bai  --attach_file_purpose bam-index           --attach_file_type bam-index
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.neg.$GENOME.bam      --attach_file_purpose neg-alignments-bam  --attach_file_type bam
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.neg.$GENOME.bam.bai  --attach_file_purpose bam-index           --attach_file_type bam-index

$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.all.${GENOME}.bw     --attach_file_purpose all-coverage-bigwig --attach_file_type bigwig
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.all.${GENOME}.starch --attach_file_purpose all-coverage-starch --attach_file_type starch
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.pos.${GENOME}.bw     --attach_file_purpose pos-coverage-bigwig --attach_file_type bigwig
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.pos.${GENOME}.starch --attach_file_purpose pos-coverage-starch --attach_file_type starch
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.neg.${GENOME}.bw     --attach_file_purpose neg-coverage-bigwig --attach_file_type bigwig
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.neg.${GENOME}.starch --attach_file_purpose neg-coverage-starch --attach_file_type starch

$ATTACH_ALIGNMENT --attach_directory ${SAMPLE_NAME}_cufflinks/ --attach_file_purpose cufflinks-directory
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}_cufflinks/transcripts.gtf         --attach_file_purpose rna-transcripts           --attach_file_type gtf
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}_cufflinks/genes.fpkm_tracking     --attach_file_purpose gene-expression-levels    --attach_file_type fpkm
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}_cufflinks/isoforms.fpkm_tracking  --attach_file_purpose isoform-expression-levels --attach_file_type fpkm
