
cd $FASTQ_DIR

FASTQ_NAME=${FLOWCELL}_${SAMPLE_NAME}

R1_NUM_FILES=$(ls ${SAMPLE_NAME}_R1_???.fastq.gz | wc -l)

if [[ "$PAIRED" == "True" ]]; then

  R2_NUM_FILES=$(ls ${SAMPLE_NAME}_R2_???.fastq.gz | wc -l)

  if [[ "$R1_NUM_FILES" -ne "$R2_NUM_FILES" ]]; then

    echo "UNEQUAL NUMBER OF FILES FOR $SAMPLE_NAME IN $FASTQ_DIR"
    exit 1
  fi
fi

UPLOAD_SCRIPT="python3 $STAMPIPES/scripts/lims/upload_data.py --attach_file_contenttype SequencingData.flowcelllane --attach_file_objectid ${FLOWCELL_LANE_ID} --attach_file_type=gzipped-fastq"

R1_FILE=${FASTQ_NAME}_R1.fastq.gz
R2_FILE=${FASTQ_NAME}_R2.fastq.gz

R1_TMP_FILE=$(mktemp)
R2_TMP_FILE=$(mktemp)

if [ -e "$R1_FILE" ]; then
  rm "$R1_FILE"
fi

if [ -e "$R2_FILE" ]; then
  rm "$R2_FILE"
fi

echo "R1: $R1_FILE"
echo "R2: $R2_FILE"

for filenum in $(seq -f "%03g" 1 $R1_NUM_FILES)
do

  echo "Adding ${filenum} to collated files"

  cat ${SAMPLE_NAME}_R1_${filenum}.fastq.gz >> $R1_TMP_FILE

  if [[ "$PAIRED" == "True" ]]; then
    R2_FILE=${FASTQ_NAME}_R2.fastq.gz
    cat ${SAMPLE_NAME}_R2_${filenum}.fastq.gz  >> $R2_TMP_FILE
  fi

done

gzip -t $R1_TMP_FILE

if [ ! -s $R1_TMP_FILE ]; then
  echo "ERROR: $R1_TMP_FILE is 0 size"
  exit 1
fi

if [ "$PAIRED" == "True" ]; then
  gzip -t $R2_TMP_FILE
  if [ ! -s $R2_TMP_FILE ]; then
    echo "ERROR: $R2_TMP_FILE is 0 size"
    exit 1
  fi
fi

rsync "$R1_TMP_FILE" "$R1_FILE"
# Files created in temp directories do not have appropriate
# permissions; make sure our collated files can be read by
# anybody
chmod 644 $R1_FILE
rm $R1_TMP_FILE
if [[ "$PAIRED" == "True" ]] ; then
  rsync "$R2_TMP_FILE" "$R2_FILE"
  chmod 644 $R2_FILE
  rm $R2_TMP_FILE
fi


$UPLOAD_SCRIPT --attach_file_purpose r1-fastq --attach_file ${R1_FILE}

if [ -e $R2_FILE ]; then
  $UPLOAD_SCRIPT --attach_file_purpose r2-fastq --attach_file ${R2_FILE}
fi
