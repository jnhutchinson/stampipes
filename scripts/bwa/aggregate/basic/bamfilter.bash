set -x -e -o pipefail

INBAM=${1:-${AGGREGATION_FOLDER}/${LIBRARY_NAME}.${GENOME}.sorted.bam}
OUTBAM=${2:-${AGGREGATION_FOLDER}/${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam}

echo Making a uniques BAM from $INBAM to $OUTBAM

cd `dirname $OUTBAM`

if [ "$UMI" = "True" ]; then
  export EXCLUDE_FLAG=1536
  echo "UMI; setting exclude flag to 1536 (QC and dups)"
else
  export EXCLUDE_FLAG=512
  echo "No UMI; setting exclude flag to 512 (QC only)"
fi

time samtools view -F $EXCLUDE_FLAG -b $INBAM > $OUTBAM
samtools index $OUTBAM
