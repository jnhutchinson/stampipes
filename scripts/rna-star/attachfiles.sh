# Uploads important files as attachments to objects on the LIMS
# Required environment names:
#  * STAMPIPES
#  * SAMPLE_NAME
#  * GENOME
#  * ALIGN_DIR
#  * OUT_DIR
#  * ALIGNMENT_ID

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

ATTACH_ALIGNMENT="python3 $UPLOAD_SCRIPT --attach_file_contenttype SequencingData.flowcelllanealignment --attach_file_objectid ${ALIGNMENT_ID}"

cd "${ALIGN_DIR}" || exit 1

$ATTACH_ALIGNMENT --attach_directory "$ALIGN_DIR" --attach_file_purpose alignment-directory
$ATTACH_ALIGNMENT --attach_file "$OUT_DIR/Aligned.toTranscriptome.out.bam" --attach_file_purpose all-alignments-bam  --attach_file_type bam

#TODO: Attach more file types
