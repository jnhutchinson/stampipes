FASTQ_NAME=${FLOWCELL}_${SAMPLE_NAME}

R1_FILE=${FASTQ_NAME}_R1.fastq.gz
R2_FILE=${FASTQ_NAME}_R2.fastq.gz

R1_FILE_FULL=${FASTQ_NAME}_R1_76.fastq.gz
R2_FILE_FULL=${FASTQ_NAME}_R2_76.fastq.gz

if [ ! -e $R1_FILE_FULL ]; then
mv $R1_FILE $R1_FILE_FULL
fi

if [ ! -e $R2_FILE_FULL ]; then
mv $R2_FILE $R2_FILE_FULL
fi

if [ ! -e $R1_FILE ]; then
zcat $R1_FILE_FULL | fastx_trimmer -i -f $START -l $END -z -Q33  > $R1_FILE
fi

if [ ! -e $R2_FILE ]; then
zcat $R2_FILE_FULL | fastx_trimmer -i -f $START -l $END -z -Q33  > $R2_FILE
fi

UPLOAD_SCRIPT="python3 $STAMPIPES/scripts/lims/upload_data.py --attach_file_contenttype SequencingData.flowcelllane --attach_file_objectid ${FLOWCELL_LANE_ID} --attach_file_type=gzipped-fastq"

$UPLOAD_SCRIPT --attach_file_purpose r1-fastq --attach_file ${R1_FILE}
$UPLOAD_SCRIPT --attach_file_purpose r2-fastq --attach_file ${R2_FILE}
$UPLOAD_SCRIPT --attach_file_purpose r1-fastq-full --attach_file ${R1_FILE_FULL}
$UPLOAD_SCRIPT --attach_file_purpose r2-fastq-full --attach_file ${R2_FILE_FULL}
