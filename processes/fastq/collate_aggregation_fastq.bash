# only support paired currently
export PAIRED=True
# Ensure that script failures have the script quit before it deletes files
set -e

cd $AGGREGATION_FOLDER

echo "Collating $AGGREGATION_FOLDER"

R1_FILE=$AGGREGATION_FOLDER/${LIBRARY_NAME}.R1.fastq.gz
R2_FILE=$AGGREGATION_FOLDER/${LIBRARY_NAME}.R2.fastq.gz

R1_NUM_FILES=`echo $R1_FILES | wc -w`
R2_NUM_FILES=`echo $R2_FILES | wc -w`

if [ "$R1_NUM_FILES" -ge "1" ]; then

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

  cat $R1_FILES > $R1_TMP_FILE
  cat $R2_FILES > $R2_TMP_FILE

  # Ensure we have valid gzipped files
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
else
  echo "No R1 files to work with"
  exit 1
fi
