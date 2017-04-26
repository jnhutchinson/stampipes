#!/bin/bash

# Dependencies
source "$MODULELOAD"
module load zlib/1.2.8
module load bedops/2.4.19
module load bedtools/2.25.0
module load bwa/0.7.12
module load jdk/1.8.0_92
module load picard/2.8.1
module load samtools/1.3
module load gcc/4.7.2
module load R/3.2.5
module load git/2.3.3
module load coreutils/8.25
module load pigz/2.3.3
module load modwt/1.0
module load hotspot2/2.0

# Load in this order specifically, currently the python3 activation
# overwrites the default "python" call, against advice
module load python/3.5.1
module load pysam/0.9.0
module load python/2.7.11

export MAX_MISMATCHES=2
export MIN_MAPPING_QUALITY=10
export MAX_INSERT_SIZE=750

JOB_BASENAME=${SAMPLE_NAME}_${FLOWCELL}_ALIGN#${ALIGNMENT_ID}

# expected output names
export FINAL_BAM=${SAMPLE_NAME}.sorted.bam
export FINAL_BAM_MARKED=${SAMPLE_NAME}.sorted.marked.bam # in case we ever want to maintain this...
export UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam
export ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
export VERSION_FILE=${SAMPLE_NAME}.versions.txt
export FASTQ_TMP=$ALIGN_DIR/fastq

cd "$ALIGN_DIR"

# check if alignment needs to be completely reset
if [[ -n "$REDO_ALIGNMENT" ]]; then
   bash "$STAMPIPES/scripts/bwa/reset_alignment.bash"
fi

# versions
bash "$STAMPIPES/scripts/versions.bash" &>"$VERSION_FILE"

# adapters
if [[ (-n "$ADAPTER_P7") && (-n "ADAPTER_P5") ]]; then
   echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" >"$ADAPTER_FILE"
fi

# trim
if [[ "$TRIM_READS_TO" -gt 0 && "$READLENGTH" -gt "$TRIM_READS_TO" ]]; then
   NEED_TRIMMING=TRUE
   READLENGTH=$TRIM_READS_TO
fi

# Indicate we have started this alignment and upload pertinent information
if [[ ! -e "$FINAL_BAM" ]]; then
   
   if [[ ! -e "$FASTQ_TMP/${SAMPLE_NAME}_R1_000.fastq.gz" ]]; then
      bash "$STAMPIPES/scripts/fastq/splitfastq.bash" "$FASTQ_TMP" "$R1_FASTQ" "$R2_FASTQ"
   fi
   
   if [[ -n "$PAIRED" ]]; then
      python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" \
         -t "$LIMS_API_TOKEN" \
         --alignment_id "$ALIGNMENT_ID" \
         --start_alignment_progress \
         --adapter_file "$ADAPTER_FILE" \
         --version_file "$VERSION_FILE"	
   else
      python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" \
         -t "$LIMS_API_TOKEN" \
  		 --start_alignment_progress \
  		 --alignment_id "$ALIGNMENT_ID" \
  		 --version_file "$VERSION_FILE"	
   fi
   
fi

# alignments
if [[ ! -e "$FINAL_BAM" ]]; then
   # individual alignments
   mkdir -p "$FASTQ_TMP"
   NUMBER_FASTQ_FILES=$(find "$FASTQ_TMP" -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l)
   FASTQ_PAIR_BAMS=""
   ALIGNMENT_JOBIDS=""
   for filenum in $(seq -f "%03g" 0 $((NUMBER_FASTQ_FILES - 1))); do
      JOBNAME=".aln${JOB_BASENAME}_${filenum}"
      BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"
      if [[ ! -e "$BAMFILE" ]]; then
         jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START ALIGNMENT: "
date

if [[ -n $PAIRED ]]; then

   fastq1=${FASTQ_TMP}/${SAMPLE_NAME}_R1_${filenum}.fastq.gz
   fastq2=${FASTQ_TMP}/${SAMPLE_NAME}_R2_${filenum}.fastq.gz

   # trim
   if [[ -n "$NEED_TRIMMING" ]] ; then
      # Define trimming function
      trim_to_length(){
         echo "Trimming \$1 to \$2 bp..."
         mv "\$1" "\$1.tmp"
         zcat \$1.tmp \
            | awk 'NR%2==0 {print substr(\$0, 1, $TRIM_TO)} NR%2!=0' \
            | gzip -c \
            > \$1
         rm \$1.tmp
      }
      trim_to_length \$fastq1 "$TRIM_READS_TO"
      trim_to_length \$fastq2 "$TRIM_READS_TO"
   fi

   # If we have UMI info, create new FASTQ files with the UMI info
   # and use those as input to the makefile
   if [[ -n "$UMI" ]]; then
      new1=\$(mktemp)
      new2=\$(mktemp)
      case "$UMI_METHOD" in
         thruplex)
            time python3 $STAMPIPES/scripts/umi/extract_umt.py \
               <(zcat \$fastq1) \
               <(zcat \$fastq2) \
               >(gzip -c -1 > \$new1) \
               >(gzip -c -1 > \$new2)
            ;;
         *)
            time python3 $STAMPIPES/scripts/umi/fastq_umi_add.py \$fastq1 \$new1
            time python3 $STAMPIPES/scripts/umi/fastq_umi_add.py \$fastq2 \$new2
            ;;
      esac
      fastq1=\$new1
      fastq2=\$new2
   fi

   # paired alignment
   make -f "$STAMPIPES/makefiles/bwa/bwa_paired_trimmed.mk" \
      "FASTQ1_FILE=\$fastq1" \
      "FASTQ2_FILE=\$fastq2" \
      "OUTBAM=$BAMFILE" \
      "TRIMSTATS=${SAMPLE_NAME}_${filenum}.trimstats.txt"

else

   # single alignment
   make -f $STAMPIPES/makefiles/bwa/bwa_single.mk \
      FASTQ_FILE=${FASTQ_TMP}/${SAMPLE_NAME}_R1_${filenum}.fastq.gz \
      OUTBAM=${BAMFILE}

fi

echo "FINISH ALIGNMENT: "
date

__SCRIPT__
)
      fi

   # Keep track of all individual alignments
   FASTQ_PAIR_BAMS="${BAMFILE} ${FASTQ_PAIR_BAMS}"
   ALIGNMENT_JOBIDS="$ALIGNMENT_JOBIDS,$jobid"

   done
   
fi

# set alignment dependencies
dependencies_align=$(echo $ALIGNMENT_JOBIDS | sed -e 's/,/,afterok:/g' | sed -e 's/^,afterok/--dependency=afterok/g')

# mark duplicates and unique file
if [[ ! -e "$FINAL_BAM" || ! -e "$UNIQUES_BAM" ]]; then

  # If we are redoing this part, then we should make sure
  # to redo all the other results as well
  export FORCE_COUNTS=1
  export REPROCESS=1
  
  JOBNAME=".pb${JOB_BASENAME}"
  export FINAL_BAM_PREFIX=${FINAL_BAM%.*} # this is for paired... because??? make file expects this??? dunno?
  finalbam_jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_align --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START BAM PROCESS: "
date

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

# mark duplicates
if [[ -n $PAIRED ]]; then
   # filter full BAM and mark duplicates
   if [[ "$UMI" == "True" ]]; then
      make -f "$STAMPIPES/makefiles/picard/dups_cigarumi.mk" SAMPLE_NAME="${SAMPLE_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
      mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
      samtools view -b -F 1536 ${FINAL_BAM} > ${UNIQUES_BAM}
   else
      make -f "$STAMPIPES/makefiles/picard/dups_cigar.mk" SAMPLE_NAME="${SAMPLE_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
      mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
      samtools view -b -F 512 ${FINAL_BAM} > ${UNIQUES_BAM}
   fi
      samtools index ${FINAL_BAM}
      samtools index ${UNIQUES_BAM}
else
   make -f $STAMPIPES/makefiles/bwa/process_unpaired_bam.mk
   make -f $STAMPIPES/makefiles/picard/dups.mk
fi

python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --flowcell_lane_id ${FLOWCELL_LANE_ID} \
   --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard

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
if [[ -n "$PAIRED" && ! -s "$SAMPLE_NAME.CollectInsertSizeMetrics.picard" ]]; then
  JOBNAME=".insert${JOB_BASENAME}"
  jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_finalbam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START INSERTSIZE: "
date

make -f "$STAMPIPES/makefiles/picard/insert_size_metrics.mk" "SAMPLE_NAME=${SAMPLE_NAME}" "BAMFILE=${UNIQUES_BAM}"

python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --flowcell_lane_id ${FLOWCELL_LANE_ID} \
   --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard \
   --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard
      
echo "END INSERTSIZE: "
date
    
__SCRIPT__
)
  PROCESSING="$PROCESSING,$jobid"
fi

# Preseq duplicate estimation
if [[ -n "$PAIRED" && ! -s "$SAMPLE_NAME.uniques.preseq.targets.txt" ]]; then
   hist=$SAMPLE_NAME.uniques.duphist.txt
   preseq=$SAMPLE_NAME.uniques.preseq.txt
   targets=$SAMPLE_NAME.uniques.preseq.targets.txt
   lorentz=$SAMPLE_NAME.uniques.dup.lorentz.txt

   JOBNAME=".ps${JOB_BASENAME}"
   jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_finalbam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START PRESEQ: "
date

# create histogram
python3 "$STAMPIPES/scripts/bam/mark_dups.py" -i "$UNIQUES_BAM" -o /dev/null --hist "$hist"

# get gini / robinhood metrics
python3 "$STAMPIPES/scripts/utility/lorentz.py" < "$hist" > "$lorentz"

# get preseq metrics
preseq lc_extrap -hist "$hist" -extrap 1.001e9 -s 1e6 -v > "$preseq"

# write preseq targets
bash "$STAMPIPES/scripts/utility/preseq_targets.sh" "$preseq" "$targets"

# uploads (condense/format these later) first one doesn't work and probably never has or at least not for awhile
#python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" -t "$LIMS_API_TOKEN" --alignment_id "$ALIGNMENT_ID" --countsfile "$lorentz"
python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" -t "$LIMS_API_TOKEN" \
   --attach_file "$preseq" \
   --attach_file_purpose preseq-estimation \
   --attach_file_type plaintext \
   --attach_file_objectid "$ALIGNMENT_ID" \
   --attach_file_contenttype SequencingData.flowcelllanealignment
python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" -t "$LIMS_API_TOKEN" --alignment_id "$ALIGNMENT_ID" --countsfile "$targets"

echo "END PRESEQ: "
date

__SCRIPT__
)
  PROCESSING="$PROCESSING,$jobid"
fi

# hotspots2
if [[ -n "$PAIRED" && ! -e "$SAMPLE_NAME.uniques.spot2.out" ]]; then
   JOBNAME=".sp2$JOB_BASENAME"
   HOTSPOT_PREFIX=$(basename "$UNIQUES_BAM" .bam)
   jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_finalbam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -e -u -o errexit

echo "Hostname: "
hostname

echo "START HOTSPOTS2: "
date

num_fragments=\$(samtools idxstats "$UNIQUES_BAM" | awk '{x+=\$3}END{print x/2}')
if [[ \$num_fragments -lt 10000000 ]]; then
   echo "Not enough fragments to get a reliable SPOT2 score, skipping"
   exit
fi
hotspot2.sh \
   -c "$BWAINDEX.chrom_sizes.bed" \
   -C "$BWAINDEX.K$READLENGTH.center_sites.n100.starch" \
   -M "$BWAINDEX.K$READLENGTH.mappable_only.bed" \
   -s 12345 \
   -F 0.1 \
   -f 0.05 \
   "$UNIQUES_BAM" \
   "$ALIGN_DIR/hotspot2"
    
cd "$ALIGN_DIR"
cp "hotspot2/$HOTSPOT_PREFIX.SPOT.txt" "$SAMPLE_NAME.uniques.spot2.out"
"$STAMPIPES/scripts/SPOT/info.sh" \
   "hotspot2/$HOTSPOT_PREFIX.hotspots.fdr0.05.starch" \
   hotspot2 \
   hotspot2/$HOTSPOT_PREFIX.SPOT.txt > "$HOTSPOT_PREFIX.hotspot2.info"

echo "END HOTSPOTS2: "
date

__SCRIPT__
)
  PROCESSING="$PROCESSING,$jobid"
fi

# tag counting
if [[ ! -e "$SAMPLE_NAME.tagcounts.txt" || -n "$FORCE_COUNTS" ]]; then
   JOBNAME=".ct${JOB_BASENAME}"
   jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_finalbam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START TAG COUNTS: "
date    

if [[ -n "$PAIRED" ]]; then
   time bash $STAMPIPES/scripts/bwa/tagcounts.bash $SAMPLE_NAME $SAMPLE_NAME.sorted.bam $SAMPLE_NAME.tagcounts.txt $R1_FASTQ $R2_FASTQ
else
   time bash $STAMPIPES/scripts/bwa/tagcounts.bash $SAMPLE_NAME $SAMPLE_NAME.sorted.bam $SAMPLE_NAME.tagcounts.txt $FASTQ_DIR
fi

# upload all data to the LIMS
python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --flowcell_lane_id ${FLOWCELL_LANE_ID} \
   --countsfile ${SAMPLE_NAME}.tagcounts.txt

echo "FINISH TAG COUNTS: "
date

__SCRIPT__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# SPOT SCORE
if [[ -n "$PAIRED" && ! -e "$SAMPLE_NAME.R1.rand.uniques.sorted.spotdups.txt" ]] || [[ ! -n "$PAIRED" && ! -e "$SAMPLE_NAME.rand.uniques.sorted.spotdups.txt" ]]; then
   JOBNAME=".sp${JOB_BASENAME}"
   jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_finalbam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START SPOT SCORE: "
date
   
if [[ -n "$PAIRED" ]]; then
   make -f $STAMPIPES/makefiles/SPOT/spot-R1-paired.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
   READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
   python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt
else
   make -f $STAMPIPES/makefiles/SPOT/spot-single.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
   READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
   python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.rand.uniques.sorted.spotdups.txt
fi

echo "FINISH SPOT SCORE: "
date
    
__SCRIPT__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# density file
if [[ ! -e "$SAMPLE_NAME.75_20.$GENOME.bw" ]]; then
   JOBNAME=".den${JOB_BASENAME}"
   jobid=$(sbatch --export=ALL -J "$JOBNAME" -o "$JOBNAME.o%A" -e "$JOBNAME.e%A" $dependencies_finalbam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START DENSITY: "
date

make -f $STAMPIPES/makefiles/densities/density.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
   READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME

echo "FINISH DENSITY: "
date

__SCRIPT__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# upload and complete
dependencies_uniquebam=$(echo $PROCESSING | sed -e 's/,/,afterok:/g' | sed -e 's/^,afterok/--dependency=afterok/g')
if [[ -n "$PROCESSING" ]]; then
  sbatch --export=ALL -J ".com$JOB_BASENAME" -o ".com$JOB_BASENAME.o%A" -e ".com$JOB_BASENAME.e%A" $dependencies_uniquebam --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=1000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START COMPLETION: "
date

if [[ -n "$PAIRED" ]]; then
   bash $STAMPIPES/scripts/bwa/checkcomplete.bash
else
   bash $STAMPIPES/scripts/bwa/unpairedcheckcomplete.bash
fi
	
python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
   -t ${LIMS_API_TOKEN} \
   -f ${FLOWCELL} \
   --alignment_id ${ALIGNMENT_ID} \
   --finish_alignment

bash $STAMPIPES/scripts/bwa/attachfiles.bash

rm -rf $ALIGN_DIR/fastq

echo "FINISH COMPLETION: "
date

__SCRIPT__
fi