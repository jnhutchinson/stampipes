source $MODULELOAD
module load bedops/2.4.19
module load jdk/1.8.0_92
module load gcc/4.7.2
module load R/3.2.5
module load picard/2.8.1
module load samtools/1.3
module load git/2.3.3
module load coreutils/8.25
module load modwt/1.0
module load hotspot2/2.0
module load bedtools/2.25.0
module load python/3.5.1
module load pysam/0.9.0

WIN=75
BINI=20
ASSAY=DNAse1 # hardcoded for now? do we need this for SPOT? isn't every assay here DNAse? (chipseq?)

export FINAL_BAM=${LIBRARY_NAME}.${GENOME}.sorted.bam
export FINAL_BAM_MARKED=${LIBRARY_NAME}.${GENOME}.sorted.marked.bam
export FINAL_UNIQUES_BAM=${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam
export TAGCOUNTS_FILE=${LIBRARY_NAME}.tagcounts.txt
export ADAPTER_COUNT_FILE=${LIBRARY_NAME}.adaptercounts.txt
export DENSITY_STARCH=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.uniques-density.bed.starch
export DENSITY_BIGWIG=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.bw
export NORM_DENSITY_BIGWIG=${LIBRARY_NAME}.${WIN}_${BINI}.normalized.${GENOME}.bw
export CUTCOUNTS_BIGWIG=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cutcounts.$READ_LENGTH.bw
export INSERT_FILE=${LIBRARY_NAME}.CollectInsertSizeMetrics.picard
export DUPS_FILE=${LIBRARY_NAME}.MarkDuplicates.picard

export HOTSPOT2_DIR=peaks
HOTSPOT_PREFIX=$(basename "$FINAL_UNIQUES_BAM" .bam)
export HOTSPOT_CALLS=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.hotspots.fdr0.05.starch
export HOTSPOT_DENSITY=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.density.bw
export HOTSPOT_SCRIPT="hotspot2.sh"
export MAPPABLE_REGIONS=${MAPPABLE_REGIONS:-$GENOME_INDEX.K${READ_LENGTH}.mappable_only.bed}
export CHROM_SIZES=${CHROM_SIZES:-$GENOME_INDEX.chrom_sizes.bed}
export CENTER_SITES=${CENTER_SITES:-$GENOME_INDEX.K${READ_LENGTH}.center_sites.n100.starch}

JOB_BASENAME=".AGG#${AGGREGATION_ID}"
MERGE_DUP_JOBNAME=${JOB_BASENAME}_merge_dup
PROCESS_BAM_JOBNAME=${JOB_BASENAME}_pb
HOTSPOT_JOBNAME=${JOB_BASENAME}_hotspot
SPOTSCORE_JOBNAME=${JOB_BASENAME}_spotscore
COUNT_JOBNAME=${JOB_BASENAME}_count
ADAPTERCOUNT_JOBNAME=${JOB_BASENAME}_adaptercount
DENSITY_JOBNAME=${JOB_BASENAME}_density
CUTCOUNTS_JOBNAME=${JOB_BASENAME}_cutcounts

cd $AGGREGATION_FOLDER
BAM_COUNT=`ls $BAM_FILES | wc -l`

if [ -n "$REDO_AGGREGATION" ]; then
        bash $STAMPIPES/scripts/bwa/aggregate/basic/reset.bash
fi

# Check out files match first
#python3 $STAMPIPES/scripts/utility/md5check.py bamfiles.txt || exit 1

PROCESSING=""

# merge bams
if [[ ! -s "$FINAL_BAM.bai" ]]; then
	merge_jobid=$(sbatch --export=ALL -J "$MERGE_DUP_JOBNAME" -o "$MERGE_DUP_JOBNAME.o%A" -e "$MERGE_DUP_JOBNAME.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: merging BAMs"
date

if [[ $BAM_COUNT -eq 1 ]]; then
	rsync ${BAM_FILES} ${FINAL_BAM}
else
	samtools merge ${FINAL_BAM} ${BAM_FILES}
fi
samtools index ${FINAL_BAM}

echo "FINISH: merging BAMs"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$merge_jobid"
fi

# set final bam dependencies
if [[ -n $merge_jobid ]]; then
   dependencies_merge=$(echo $merge_jobid | sed -e 's/^/--dependency=afterok:/g')
fi
	
# process BAM file
if [[ ! -s "$FINAL_UNIQUES_BAM.bai" ]]; then	
	pb_jobid=$(sbatch --export=ALL -J "$MERGE_DUP_JOBNAME" -o "$MERGE_DUP_JOBNAME.o%A" -e "$MERGE_DUP_JOBNAME.e%A" $dependencies_merge --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: processing BAM"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

if [[ "$UMI" == "True" && -n "$PAIRED" ]]; then
	make -f "$STAMPIPES/makefiles/picard/dups_cigarumi.mk" SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
	mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
	samtools view -b -F 1536 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
elif [[ -n "$PAIRED" ]]; then
	make -f "$STAMPIPES/makefiles/picard/dups_cigar.mk" SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
	mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
	samtools view -b -F 512 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
else
        make -f $STAMPIPES/makefiles/picard/dups.mk SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM=${FINAL_BAM_MARKED}
        mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
        samtools view -b -F 512 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
fi

samtools index ${FINAL_BAM}
samtools index ${FINAL_UNIQUES_BAM}

# calculate insert sizes
if [[ ! -s "$INSERT_FILE" && -n "$PAIRED" ]]; then
	make -f "$STAMPIPES/makefiles/picard/insert_size_metrics.mk" "SAMPLE_NAME=${LIBRARY_NAME}" "BAMFILE=${FINAL_UNIQUES_BAM}" "INSERTMETRICS=${INSERT_FILE}"
fi

rm -rf "\$TMPDIR"

echo "FINISH: processing BAM"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$pb_jobid"
fi

# set final bam dependencies
if [[ -n $pb_jobid ]]; then
   dependencies_pb=$(echo $pb_jobid | sed -e 's/^/--dependency=afterok:/g')
fi

# Run Hotspot2
if [[ ! -s "$HOTSPOT_CALLS" || ! -s "$HOTSPOT_DENSITY" ]] ; then
	HOTSPOT_SPOT=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.SPOT.txt
	jobid=$(sbatch --export=ALL -J "$HOTSPOT_JOBNAME" -o "$HOTSPOT_JOBNAME.o%A" -e "$HOTSPOT_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

echo "Hostname: "
hostname

echo "START: hotspot2"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

"$HOTSPOT_SCRIPT"  -F 0.5 -s 12345 -M "$MAPPABLE_REGIONS" -c "$CHROM_SIZES" -C "$CENTER_SITES" "$FINAL_UNIQUES_BAM"  "$HOTSPOT2_DIR"
"$STAMPIPES/scripts/SPOT/info.sh" "$HOTSPOT_CALLS" hotspot2 \$(cat $HOTSPOT_SPOT) > "$HOTSPOT_PREFIX.hotspot2.info"

echo "FINISH: hotspot2"
date

rm -rf "\$TMPDIR"

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# SPOT score
if [[ -n "$PAIRED" && ! -e "$LIBRARY_NAME.$GENOME.R1.rand.uniques.sorted.spotdups.txt" ]] || [[ ! -n "$PAIRED" && ! -e "$LIBRARY_NAME.$GENOME.rand.uniques.sorted.spotdups.txt" ]]; then
	jobid=$(sbatch --export=ALL -J "$SPOTSCORE_JOBNAME" -o "$SPOTSCORE_JOBNAME.o%A" -e "$SPOTSCORE_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: HOTSPOT1 SPOT SCORE"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

if [[ -n "$PAIRED" ]]; then
	make -f $STAMPIPES/makefiles/SPOT/spot-R1-paired.mk BWAINDEX=$GENOME_INDEX ASSAY=$ASSAY GENOME=$GENOME READLENGTH=$READ_LENGTH SAMPLE_NAME="$LIBRARY_NAME.$GENOME"
else
	make -f $STAMPIPES/makefiles/SPOT/spot-single.mk BWAINDEX=$GENOME_INDEX ASSAY=$ASSAY GENOME=$GENOME READLENGTH=$READ_LENGTH SAMPLE_NAME="$LIBRARY_NAME.$GENOME"
fi

rm -rf "\$TMPDIR"

echo "FINISH: HOTSPOT1 SPOT SCORE"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# bam counts
if [ ! -e $TAGCOUNTS_FILE ]; then
	jobid=$(sbatch --export=ALL -J "$COUNT_JOBNAME" -o "$COUNT_JOBNAME.o%A" -e "$COUNT_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: tag counts"
date

python3 $STAMPIPES/scripts/bwa/bamcounts.py $FINAL_BAM $TAGCOUNTS_FILE

echo "FINISH: tag counts"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# adapter counts
if [[ ! -s "$ADAPTER_COUNT_FILE" ]]; then
	jobid=$(sbatch --export=ALL -J "$ADAPTERCOUNT_JOBNAME" -o "$ADAPTERCOUNT_JOBNAME.o%A" -e "$ADAPTERCOUNT_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: adapter counts"
date

adaptercount=\$(bash "$STAMPIPES/scripts/bam/count_adapters.sh" "$FINAL_BAM")
if [ -n \$adapter_count ]; then
	echo -e "adapter\t\$adaptercount" > "$ADAPTER_COUNT_FILE"
fi

echo "FINISH: adapter counts"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# density tracks
if [[ ! -s "$DENSITY_BIGWIG" || ! -s "$NORM_DENSITY_BIGWIG" ]]; then
	jobid=$(sbatch --export=ALL -J "$DENSITY_JOBNAME" -o "$DENSITY_JOBNAME.o%A" -e "$DENSITY_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: density"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

make -f $STAMPIPES/makefiles/densities/density.mk FAI=${GENOME_INDEX}.fai SAMPLE_NAME=${LIBRARY_NAME} GENOME=${GENOME} \
	BAMFILE=${FINAL_UNIQUES_BAM} STARCH_OUT=${DENSITY_STARCH} BIGWIG_OUT=${DENSITY_BIGWIG}
make -f "$STAMPIPES/makefiles/densities/normalize-density.mk" BAMFILE=${FINAL_UNIQUES_BAM} SAMPLE_NAME=${LIBRARY_NAME} FAI=${GENOME_INDEX}.fai

rm -rf "\$TMPDIR"

echo "FINISH: density"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# cutcounts
if [ ! -e "$CUTCOUNTS_BIGWIG" ]; then
	jobid=$(sbatch --export=ALL -J "$CUTCOUNTS_JOBNAME" -o "$CUTCOUNTS_JOBNAME.o%A" -e "$CUTCOUNTS_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: cutcounts"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

bash $STAMPIPES/scripts/bwa/aggregate/basic/cutcounts.bash

rm -rf "\$TMPDIR"

echo "FINISH: cutcounts"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# get complete dependencies
complete_dependencies=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterok/g')

# upload data
if [[ -n "${PROCESSING}" ]]; then
	UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py
	jobid=$(sbatch --export=ALL -J "${JOB_BASENAME}_complete" -o "${JOB_BASENAME}_complete.o%A" -e "${JOB_BASENAME}_complete.e%A" $complete_dependencies --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=1000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: "
date

bash "$STAMPIPES/scripts/bwa/aggregate/basic/checkcomplete.bash"
bash "$STAMPIPES/scripts/bwa/aggregate/basic/attachfiles.bash"
bash "$STAMPIPES/scripts/bwa/aggregate/basic/uploadcounts.bash"

echo "FINISH: "
date

__SCRIPT__
)
fi