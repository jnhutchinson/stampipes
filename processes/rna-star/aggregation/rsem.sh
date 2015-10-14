source $MODULELOAD
module load samtools/1.2
module load gcc/4.7.2     # R dependent on this
module load R/3.1.0

module load RSEM/1.2.22

export REFDIR="$(dirname $GENOME_INDEX)"
export RSEMrefDir="$REFDIR/RSEMgenome-hg19-g19-combined/"

export TARGET_BAM=Aligned.toTranscriptome.out.bam

numbam=$(wc -w <<< $BAM_FILES)
# Temporary
if [ ! -s "$TARGET_BAM" ] ; then
  if [ $numbam -eq 1 ] ; then
    cp "$BAM_FILES" "$TARGET_BAM"
  else
    samtools merge -n -f "$TARGET_BAM" $BAM_FILES
  fi
fi

if [ ! -s "Quant.genes.results" ] ; then

  qsub -cwd -V -N ".AGG#$AGGREGATION_ID" -pe threads 2-4 -S /bin/bash <<'__RSEM__'
    nThreadsRSEM=$((NSLOTS * 2))

    SORTED_BAM=$TARGET_BAM

    RSEM=rsem-calculate-expression

    # RSEM parameters: common
    RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345 --temporary-folder $TMPDIR"

    # RSEM parameters: run-time, number of threads and RAM in MB
    RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

    # RSEM parameters: data type dependent
    #OPTION: stranded paired end
    RSEMparType="--paired-end --forward-prob 1.0"

    RSEMrefDir=$($STAMPIPES/scripts/cache.sh $RSEMrefDir)

    ###### RSEM command
    (
    echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType $SORTED_BAM $RSEMrefDir/RSEMref Quant
         $RSEM $RSEMparCommon $RSEMparRun $RSEMparType $SORTED_BAM $RSEMrefDir/RSEMref Quant
    ) >& Log.rsem

    echo rsem-plot-model Quant Quant.pdf
    rsem-plot-model Quant Quant.pdf
__RSEM__

fi
