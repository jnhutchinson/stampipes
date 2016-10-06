source $MODULELOAD
module load samtools/1.3
module load gcc/4.7.2     # R dependent on this
module load R/3.2.5

module load STAR/2.4.2a   # Just for densities
module load bedops/2.4.19

module load RSEM/1.2.30

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
  samtools index "$GENOME_BAM"
fi

if [ ! -s "Signal.UniqueMultiple.str+.starch" ] ; then
  qsub -cwd -V -N ".AGG#${AGGREGATION_ID}.den" -S /bin/bash <<'__DEN__'

    # Write starch and bigwig to .tmp files
    function convertBedGraph(){
      in="$1"
      base="$2"
      chrom="$in.onlyChr.bg"
      grep '^chr' "$in" > $chrom
      bedGraphToBigWig "$chrom" chrNL.txt "$base.bw.tmp"
      starch "$chrom" > "$base.starch.tmp"
    }

    set -x -e -o pipefail

    mkdir -p $TMPDIR/Signal

    echo STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR

    grep '^chr' $STARrefDir/chrNameLength.txt > chrNL.txt

    convertBedGraph $TMPDIR/Signal/Signal.Unique.str1.out.bg         Signal.Unique.str-
    convertBedGraph $TMPDIR/Signal/Signal.Unique.str2.out.bg         Signal.Unique.str+
    convertBedGraph $TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str-
    convertBedGraph $TMPDIR/Signal/Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str+

    for i in Signal*.tmp ; do
      mv $i ${i/.tmp/}
    done

__DEN__

fi

if [ ! -s "Quant.genes.results" ] ; then

  qsub -cwd -V -N ".AGG#${AGGREGATION_ID}.rsem" -pe threads 4 -S /bin/bash <<'__RSEM__'
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
