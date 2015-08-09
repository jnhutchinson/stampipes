# Environment variables expected:
#
# * FASTQ_TMP: where to put the fastq
# * SAMPLE_NAME: the "base name" of the FASTQ files
# * R1_FASTQ: The fastq file for Read 1
# * R2_FASTQ: The fastq file for Read 2

R1_FASTQ=$1
R2_FASTQ=$2
FASTQ_TMP=$3
CHUNK_SIZE=16000000

LINE_CHUNK=$(($CHUNK_SIZE * 4))

if [[ -e $FASTQ_TMP ]]; then
  rm -rf $FASTQ_TMP/${SAMPLE_NAME}_R?_???.fastq.gz
fi

mkdir -p $FASTQ_TMP

echo "Splitting FASTQ files into $CHUNK_SIZE reads using $TMPDIR"
date

echo "Splitting $R1_FASTQ"
zcat $R1_FASTQ | split -l $LINE_CHUNK -d -a 3 - "$TMPDIR/${SAMPLE_NAME}_R1_"
date
echo "Splitting $R2_FASTQ"
zcat $R2_FASTQ | split -l $LINE_CHUNK -d -a 3 - "$TMPDIR/${SAMPLE_NAME}_R2_"
date

SPLIT_COUNT=`ls $TMPDIR/${SAMPLE_NAME}_R1_??? | wc -l`

if [ "$SPLIT_COUNT" -eq 1 ]; then
    # We only have one file split
    # Don't bother compressing again, just symlink the existing full fastq files
    echo "Files do not exceed $CHUNK_SIZE; symlinking original files"
    ln -s $R1_FASTQ $FASTQ_TMP/${SAMPLE_NAME}_R1_000.fastq.gz
    ln -s $R2_FASTQ $FASTQ_TMP/${SAMPLE_NAME}_R2_000.fastq.gz
    exit
fi

for RAW_FILE in `find ${TMPDIR} -name "${SAMPLE_NAME}_R?_???"`; do
    mv $RAW_FILE ${RAW_FILE}.fastq
    FASTQ_FILENAME=`basename $RAW_FILE`
    pigz --fast -p 2 -c ${RAW_FILE}.fastq > $FASTQ_TMP/$FASTQ_FILENAME.fastq.gz
    date
done
