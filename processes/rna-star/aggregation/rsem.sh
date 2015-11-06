source $MODULELOAD
module load samtools/1.2
module load gcc/4.7.2     # R dependent on this
module load R/3.1.0

module load star/2.4.2a   # Just for densities
module load RSEM/1.2.22

export REFDIR="$(dirname $GENOME_INDEX)"
export STARrefDir="$REFDIR/STARgenome-hg19-g19-combined/"
export RSEMrefDir="$REFDIR/RSEMgenome-hg19-g19-combined/"

export TARGET_BAM=Aligned.toTranscriptome.out.bam
export GENOME_BAM=Aligned.toGenome.out.bam

numbam=$(wc -w <<< $BAM_FILES)
# Temporary
if [ ! -s "$TARGET_BAM" ] ; then
  if [ $numbam -eq 1 ] ; then
    cp "$BAM_FILES" "$TARGET_BAM"
  else
    samtools merge -n -f "$TARGET_BAM" $BAM_FILES
  fi
fi

if [ ! -s "$GENOME_BAM" ] ; then
  GENOME_BAM_FILES=$(sed 's/toTranscriptome/sortedByCoord/g' <<< "$BAM_FILES")
  $STAMPIPES/scripts/tophat/merge_or_copy_bam.sh "$GENOME_BAM" $GENOME_BAM_FILES
fi

if [ ! -s "Signal.Unique.strand+.bw" ] ; then
  qsub -cwd -V -N ".AGG#${AGGREGATION_ID}.den" -S /bin/bash <<'__DEN__'
    set -x -e -o pipefail

    mkdir -p $TMPDIR/Signal

    echo STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR

    grep '^chr' $STARrefDir/chrNameLength.txt > chrNL.txt

    tree -s $TMPDIR #debug
    for i in $(ls $TMPDIR/Signal/Signal*bg); do
      grep '^chr' $i > $i.onlyChr.bg
    done

    bedGraphToBigWig $TMPDIR/Signal/Signal.Unique.str1.out.bg.onlyChr.bg         chrNL.txt Signal.Unique.str-.bw.tmp
    bedGraphToBigWig $TMPDIR/Signal/Signal.Unique.str2.out.bg.onlyChr.bg         chrNL.txt Signal.Unique.str+.bw.tmp
    bedGraphToBigWig $TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg.onlyChr.bg chrNL.txt Signal.UniqueMultiple.str-.bw.tmp
    bedGraphToBigWig $TMPDIR/Signal/Signal.UniqueMultiple.str2.out.bg.onlyChr.bg chrNL.txt Signal.UniqueMultiple.str+.bw.tmp

    for i in Signal*bw.tmp ; do
      mv $i ${i/.tmp/}
    done

__DEN__

fi

if [ ! -s "Quant.genes.results" ] ; then

  qsub -cwd -V -N ".AGG#${AGGREGATION_ID}.rsem" -pe threads 2-4 -S /bin/bash <<'__RSEM__'
    set -x

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
