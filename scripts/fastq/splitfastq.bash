# Environment variables expected:
#
# * FASTQ_TMP: where to put the fastq
# * SAMPLE_NAME: the "base name" of the FASTQ files
# * R1_FASTQ: The fastq file for Read 1
# * R2_FASTQ: The fastq file for Read 2

CHUNK_SIZE=16000000

LINE_CHUNK=$(($CHUNK_SIZE * 4))

if [[ -e $FASTQ_TMP ]]; then
  rm -rf $FASTQ_TMP/${SAMPLE_NAME}_R?_???.fastq.gz
fi

mkdir -p $FASTQ_TMP

echo "Splitting FASTQ files into $CHUNK_SIZE reads"

echo "Splitting $R1_FASTQ"
zcat $R1_FASTQ | split -l $LINE_CHUNK -d -a 3 - "$FASTQ_TMP/${SAMPLE_NAME}_R1_"
echo "Splitting $R2_FASTQ"
zcat $R2_FASTQ | split -l $LINE_CHUNK -d -a 3 - "$FASTQ_TMP/${SAMPLE_NAME}_R2_"

for RAW_FILE in `find ${FASTQ_TMP} -name "${SAMPLE_NAME}_R?_???"`; do
    mv $RAW_FILE ${RAW_FILE}.fastq
    gzip ${RAW_FILE}.fastq
done
