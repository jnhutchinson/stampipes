# Uploads important files as attachments to an aggregation object on the LIMS
# Required environment names:
#  * STAMPIPES
#  * LIBRARY_NAME
#  * GENOME
#  * AGGREGATION_DIR
#  * AGGREGATION_ID

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

ATTACH_AGGREGATION="python3 $UPLOAD_SCRIPT --attach_file_contenttype AggregationData.aggregation --attach_file_objectid ${AGGREGATION_ID}"

$ATTACH_AGGREGATION --attach_directory `pwd` --attach_file_purpose aggregation-directory 
$ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.${GENOME}.sorted.bam --attach_file_purpose all-alignments-bam --attach_file_type bam
$ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.${GENOME}.sorted.bam.bai --attach_file_purpose bam-index --attach_file_type bam-index

# Densities
$ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.75_20.${GENOME}.bw --attach_file_purpose density-bigwig-windowed --attach_file_type bigwig
$ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.75_20.${GENOME}.uniques-density.bed.starch --attach_file_purpose density-bed-starch-windowed --attach_file_type starch
$ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.75_20.normalized.${GENOME}.bw --attach_file_purpose normalized-density-bigwig-windowed --attach_file_type bigwig
$ATTACH_AGGREGATION --attach_file ${LIBRARY_NAME}.75_20.normalized.${GENOME}.uniques-density.bed.starch --attach_file_purpose normalized-density-bed-starch --attach_file_type starch
