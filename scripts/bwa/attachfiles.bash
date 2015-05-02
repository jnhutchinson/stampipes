# Uploads important files as attachments to objects on the LIMS
# Required environment names:
#  * STAMPIPES
#  * SAMPLE_NAME
#  * GENOME
#  * ALIGN_DIR
#  * ALIGNMENT_ID
#  * READLENGTH

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

ATTACH_ALIGNMENT="python $UPLOAD_SCRIPT --attach_file_contenttype SequencingData.flowcelllanealignment --attach_file_objectid ${ALIGNMENT_ID}"

$ATTACH_ALIGNMENT --attach_directory ${ALIGN_DIR} --attach_file_purpose alignment-directory 
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.sorted.bam --attach_file_purpose all-alignments-bam --attach_file_type bam
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.sorted.bam.bai --attach_file_purpose bam-index --attach_file_type bam-index
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.uniques.sorted.bam --attach_file_purpose filtered-alignments-bam --attach_file_type bam
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.uniques.sorted.bam.bai --attach_file_purpose bam-index --attach_file_type bam-index
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.75_20.${GENOME}.bw --attach_file_purpose density-bigwig-windowed --attach_file_type bigwig
$ATTACH_ALIGNMENT --attach_file ${SAMPLE_NAME}.75_20.uniques-density.${READLENGTH}.${GENOME}.bed.starch --attach_file_purpose density-bed-starch-windowed --attach_file_type starch

