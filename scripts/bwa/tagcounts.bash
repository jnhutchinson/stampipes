SAMPLE_NAME=$1
INBAM=$2
OUTPUT=$3
R1_FASTQ=$4
R2_FASTQ=$5

# Delete output if it already exists; do not want to append to existing output
if [ -e $OUTPUT ]; then
    rm $OUTPUT
fi

if [ ! -e $INBAM ]; then
  echo "Could not find $INBAM"
  exit 1
fi

echo "Calculating total/pf/qc counts"
if [[ "$PAIRED" == "True" ]]; then
    AWK_PAIRED="-v paired=1";
else
    AWK_PAIRED="-v paired=0";
fi
if [ -e $R1_FASTQ ]; then
  echo $AWK_PAIRED
  zcat $R1_FASTQ | awk $AWK_PAIRED -f $STAMPIPES/awk/illumina_fastq_count.awk >> $OUTPUT
fi

echo "Calculate tags trimmed"
# count is in read pairs, double for accuracy in LIMS display % wise
TRIMMED=`find . -maxdepth 1 -name "$SAMPLE_NAME*trimstats.txt" | xargs awk 'BEGIN { total=0 } { match($0, /Total read-pairs trimmed: ([0-9]+)/, a); total=total+a[1];} END{ print total * 2 }'`
echo -e "adapter-trimmed\t$TRIMMED" >> $OUTPUT

echo "Creating bam counts"
if [[ "$PAIRED" = "True" ]]; then
  python3 $STAMPIPES/scripts/bwa/bamcounts.py $INBAM $OUTPUT
else
  python3 $STAMPIPES/scripts/bwa/bamcounts.py --unpaired $INBAM $OUTPUT
fi
