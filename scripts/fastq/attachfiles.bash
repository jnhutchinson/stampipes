# Uploads important files as attachments to objects on the LIMS
# Required environment names:
#  * STAMPIPES
#  * SAMPLE_NAME
#  * GENOME
#  * FASTQ_DIR
#  * FLOWCELL_LANE_ID

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

ATTACH_LANE="python3 $UPLOAD_SCRIPT --attach_file_contenttype SequencingData.flowcelllane --attach_file_object ${FLOWCELL_LANE_ID}"
ATTACH_FASTQ="$ATTACH_LANE --attach_file_purpose raw-reads --attach_file_type gzipped-fastq"

$ATTACH_LANE --attach_directory ${FASTQ_DIR} --attach_file_purpose fastq-directory

$ATTACH_LANE --attach_file ${FASTQ_DIR}/${SAMPLE_NAME}_R1_fastqc.zip --attach_file_type zip --attach_file_purpose fastqc-results-zip

if [ "$PAIRED" = "True" ]; then
    $ATTACH_LANE --attach_file ${FASTQ_DIR}/${SAMPLE_NAME}_R2_fastqc.zip --attach_file_type zip --attach_file_purpose fastqc-results-zip
fi

for FASTQ_FILE in `find $FASTQ_DIR -name "$SAMPLE_NAME*.fastq.gz"`; do
    $ATTACH_FASTQ --attach_file $FASTQ_FILE
done
