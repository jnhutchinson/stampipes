#!/bin/bash

# module information
source "$MODULELOAD"
module load bwa/0.7.12
module load samtools/1.3
module load jdk/1.8.0_92
module load picard/2.8.1

# job suffix
JOB_BASENAME=${SAMPLE_NAME}_${FLOWCELL}_ALIGN#${ALIGNMENT_ID}

# mapping parameters
export THREAD_NUM=3
export MAX_MISMATCHES=2
export MIN_MAPPING_QUALITY=10
export MAX_INSERT_SIZE=750
export BWA_ALN_PARAM="-t 3 -Y -l 32 -n 0.04"
export BWA_SAMPE_PARAM="-n 1000000 -a 750"

# expected output names
export ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
export VERSION_FILE=${SAMPLE_NAME}.versions.txt
export FASTQ_TMP=$ALIGN_DIR/fastq
export FINAL_BAM=${SAMPLE_NAME}.sorted.bam
export FINAL_UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam
export INSERT_FILE=${SAMPLE_NAME}.CollectInsertSizeMetrics.picard
export INSERT_HIST=${SAMPLE_NAME}.CollectInsertSizeMetrics.picard.pdf
export NUCLEAR_CHR=${NUCLEAR_CHR:-$BWAINDEX.nuclear.txt}

###

cd "$ALIGN_DIR"

# versions
bash "$STAMPIPES/scripts/versions.bash" &>"$VERSION_FILE"

# adapters
if [[ (-n "$ADAPTER_P7") && (-n "ADAPTER_P5") ]]; then
   echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" >"$ADAPTER_FILE"
fi

# alignment
if [[ ! -e "$FINAL_BAM" ]]; then

    # split FASTQ
   if [[ ! -e "$FASTQ_TMP/${SAMPLE_NAME}_R1_000.fastq.gz" ]]; then
      bash "$STAMPIPES/scripts/fastq/splitfastq.bash" "$FASTQ_TMP" "$R1_FASTQ" "$R2_FASTQ"
   fi

   # individual alignments
   mkdir -p "$FASTQ_TMP"
   NUMBER_FASTQ_FILES=$(find "$FASTQ_TMP" -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l)
   FASTQ_PAIR_BAMS=""
   ALIGNMENT_JOBIDS=""
   
   # go through each split FASTQ and align individually
   for filenum in $(seq -f "%03g" 0 $((NUMBER_FASTQ_FILES - 1))); do
      JOBNAME=".aln${JOB_BASENAME}_${filenum}"
      BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"
      filenum=$filenum
      if [[ ! -e "$BAMFILE" ]]; then
         jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START ALIGNMENT: "
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

fastq1=${FASTQ_TMP}/${SAMPLE_NAME}_R1_${filenum}.fastq.gz
fastq2=${FASTQ_TMP}/${SAMPLE_NAME}_R2_${filenum}.fastq.gz

time trim-adapters-illumina -f $ADAPTER_FILE -1 P5 -2 P7 --threads=$THREAD_NUM $fastq1 $fastq2 \$TMPDIR/trimmed.R1.fastq.gz \$TMPDIR/trimmed.R2.fastq.gz \
&> ${SAMPLE_NAME}_${filenum}.trimstats.txt \
&& echo trimmed ${SAMPLE_NAME} >&2
time bwa aln $BWA_ALN_PARAM $BWAINDEX \$TMPDIR/trimmed.R1.fastq.gz > \$TMPDIR/R1.sai && echo made \$TMPDIR/R1.sai >&2
time bwa aln $BWA_ALN_PARAM $BWAINDEX \$TMPDIR/trimmed.R2.fastq.gz > \$TMPDIR/R2.sai && echo made \$TMPDIR/R2.sai >&2
time bwa sampe $BWA_SAMPE_PARAM $BWAINDEX \$TMPDIR/R1.sai \$TMPDIR/R2.sai \$TMPDIR/trimmed.R1.fastq.gz \$TMPDIR/trimmed.R2.fastq.gz > \$TMPDIR/align.sam && echo made \$TMPDIR/align.sam >&2

time samtools view -bS -t ${BWAINDEX}.fai \$TMPDIR/align.sam > \$TMPDIR/align.unsorted.bam && echo made \$TMPDIR/align.unsorted.bam >&2
time python3 $STAMPIPES/scripts/bwa/filter_reads.py \$TMPDIR/align.unsorted.bam \$TMPDIR/align.filtered.bam $NUCLEAR_CHR
time samtools sort -l 0 -m 1G -@ 8 \$TMPDIR/align.filtered.bam > \$TMPDIR/align.sorted.bam && echo made \$TMPDIR/align.sorted.bam >&2
rsync \$TMPDIR/align.sorted.bam $BAMFILE

rm -rf "\$TMPDIR"

echo "FINISH ALIGNMENT: "
date

__SCRIPT__
)
      fi

   # Keep track of all individual alignments
   FASTQ_PAIR_BAMS="${BAMFILE} ${FASTQ_PAIR_BAMS}"
   if [[ -n $jobid ]]; then
       ALIGNMENT_JOBIDS="$ALIGNMENT_JOBIDS,$jobid"
   fi

done

fi

# set alignment dependencies
dependencies_align=$(echo $ALIGNMENT_JOBIDS | sed -e 's/,/,afterok:/g' | sed -e 's/^,afterok/--dependency=afterok/g')

# mark duplicates and unique file
if [[ ! -e "$FINAL_BAM" || ! -e "$FINAL_UNIQUES_BAM" ]]; then

  JOBNAME=".pb${JOB_BASENAME}"
  finalbam_jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_align --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START BAM PROCESS: "
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

# merge BAMs
if [[ ! -e ${FINAL_BAM} ]]; then
   if [[ "$NUMBER_FASTQ_FILES" -eq "1" ]]; then
      mv ${SAMPLE_NAME}_000.sorted.bam ${FINAL_BAM}
   else
      samtools merge ${FINAL_BAM} ${FASTQ_PAIR_BAMS}
   fi
   samtools index ${FINAL_BAM}
else
   echo "$FINAL_BAM exists already"
fi
if [[ "$NUMBER_FASTQ_FILES" -gt "1" ]]
then
   rm -f $FASTQ_PAIR_BAMS
fi

rm -rf "\$TMPDIR"

echo "FINISH BAM PROCESS: "
date

__SCRIPT__
)
  PROCESSING="$PROCESSING,$finalbam_jobid"
fi

# set final bam dependencies
if [[ -n $finalbam_jobid ]]; then
   dependencies_finalbam=$(echo $finalbam_jobid | sed -e 's/^/--dependency=afterok:/g')
fi

# insertsize
if [[ -n "$PAIRED" && ! -s "$INSERT_FILE" ]]; then
  JOBNAME=".insert${JOB_BASENAME}"
  jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_finalbam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START INSERTSIZE: "
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

samtools idxstats ${SAMPLE_NAME}.uniques.sorted.bam | cut -f 1 | grep -v chrM | grep -v chrC | xargs samtools view -b ${SAMPLE_NAME}.uniques.sorted.bam > ${TMPDIR}/${SAMPLE_NAME}.nuclear.bam

picard CollectInsertSizeMetrics INPUT=${TMPDIR}/${SAMPLE_NAME}.nuclear.bam OUTPUT=$INSERT_FILE \
    HISTOGRAM_FILE=$INSERT_HIST \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true && echo Picard stats >&2

python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --flowcell_lane_id ${FLOWCELL_LANE_ID} \
   --insertsfile $INSERT_FILE

rm -rf "\$TMPDIR"

echo "END INSERTSIZE: "
date

__SCRIPT__
)
  PROCESSING="$PROCESSING,$jobid"
fi

# upload and complete
dependencies_uniquebam=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterok/g')
if [[ -n "$PROCESSING" ]]; then
  upload_id=$(sbatch --export=ALL -J ".com$JOB_BASENAME" -o ".com$JOB_BASENAME.o%A" -e ".com$JOB_BASENAME.e%A" $dependencies_uniquebam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=1000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

sleep 10

__SCRIPT__
)
fi

echo "$upload_id" > last_complete_job_id.txt