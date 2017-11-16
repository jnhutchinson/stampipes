source $MODULELOAD
module load samtools/1.3
module load gcc/4.7.2     # R dependent on this
module load R/3.2.5
module load STAR/2.4.2a   # Just for densities
module load bedops/2.4.19
module load subread/1.5.1 # for featureCounts
module load cufflinks/2.2.1 # for cuffLinks
module load anaquin/2.0.1
module load kallisto/0.43.1
module load htslib/1.6.0
module load jdk/1.8.0_92
module load picard/2.8.1
source "$PYTHON3_ACTIVATE"
module load python/2.7.11

export REFDIR="$(dirname $GENOME_INDEX)"
export STARrefDir="$REFDIR/${STAR_DIR}"
export TARGET_BAM=Aligned.toTranscriptome.out.bam
export GENOME_BAM=Aligned.toGenome.out.bam

if [ -n "$REDO_AGGREGATION" ]; then
    bash $STAMPIPES/scripts/rna-star/aggregate/reset.bash
fi

TRIMS_R1=trims.R1.fastq.gz
TRIMS_R2=trims.R2.fastq.gz
# create merged fastqs
if [ ! -s "$TRIMS_R1" ] ; then
    cat $TRIMMED_R1 > $TRIMS_R1
    cat $TRIMMED_R2 > $TRIMS_R2
fi

# create proper merged BAM
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
# anaquin call on merged BAM
anaquin RnaAlign -rgtf $SEQUINS_REF -usequin $GENOME_BAM -o anaquin_star
bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaalign_stats.bash anaquin_star/RnaAlign_summary.stats anaquin_star/RnaAlign_summary.stats.info

density_job=.AG${AGGREGATION_ID}.star_den
cufflinks_job=.AG${AGGREGATION_ID}.star_cuff
kallisto_job=.AGG${AGGREGATION_ID}.star_kallisto
complete_job=.AGG#${AGGREGATION_ID}.complete
fcounts_job=.AGG${AGGREGATION_ID}.star_fcounts
picard_job=.AGG${AGGREGATION_ID}.picard
tagcounts_job=.AGG${AGGREGATION_ID}.tagcounts
adaptercounts_job=.AGG${AGGREGATION_ID}.adapter

# density information, convoluted, can clean up so we skip a lot of these steps
if [ ! -s "Signal.UniqueMultiple.str+.starch" ] ; then
    jobid=$(sbatch --export=ALL -J "$density_job" -o "$density_job.o%A" -e "$density_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__DEN__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START DENSITY: "
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

    # Write starch and bigwig to .tmp files
    function convertBedGraph(){
      in="\$1"
      base="\$2"
      chrom="\$in.onlyChr.bg"
      grep '^chr' "\$in" | sort -k1,1 -k2,2n > \$chrom
      bedGraphToBigWig "\$chrom" chrNL.txt "\$base.bw.tmp"
      starch "\$chrom" > "\$base.starch.tmp"
    }

    mkdir -p \$TMPDIR/Signal

    echo STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix \$TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir \$TMPDIR/STAR
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Unstranded --outFileNamePrefix \$TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir \$TMPDIR/STAR
    mv \$TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg \$TMPDIR/Signal/Signal.UniqueMultiple.unstranded.out.bg
    mv \$TMPDIR/Signal/Signal.Unique.str1.out.bg \$TMPDIR/Signal/Signal.Unique.unstranded.out.bg
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix \$TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir \$TMPDIR/STAR

    grep '^chr' $STARrefDir/chrNameLength.txt > chrNL.txt

    convertBedGraph \$TMPDIR/Signal/Signal.Unique.str1.out.bg         Signal.Unique.str-
    convertBedGraph \$TMPDIR/Signal/Signal.Unique.str2.out.bg         Signal.Unique.str+
    convertBedGraph \$TMPDIR/Signal/Signal.Unique.unstranded.out.bg         Signal.Unique.both
    convertBedGraph \$TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str-
    convertBedGraph \$TMPDIR/Signal/Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str+
    convertBedGraph \$TMPDIR/Signal/Signal.UniqueMultiple.unstranded.out.bg Signal.UniqueMultiple.both

    for i in Signal*.tmp ; do
      mv \$i \${i/.tmp/}
    done

    unstarch Signal.Unique.both.starch | bgzip > Signal.Unique.both.starch.bgz
    tabix -p bed Signal.Unique.both.starch.bgz

rm -rf "\$TMPDIR"

echo "FINISH DENSITY: "
date

__DEN__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# cufflinks
if [ ! -s "anaquin_cufflinks/RnaExpression_summary.stats" ] ; then
    jobid=$(sbatch --export=ALL -J "$cufflinks_job" -o "$cufflinks_job.o%A" -e "$cufflinks_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__CUFF__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START CUFFLINKS: "
date

    CUFF=cufflinks
    CUFF_COMMON="--no-update-check --library-type fr-firststrand"

    \$CUFF \$CUFF_COMMON --GTF $ANNOTATION $GENOME_BAM

    # delete duplicate rows and sort into identical orders across all samples
    Rscript $STAMPIPES/scripts/rna-star/aggregate/dedupe_sort_cuffout.Rscript genes.fpkm_tracking
    Rscript $STAMPIPES/scripts/rna-star/aggregate/dedupe_sort_cuffout.Rscript isoforms.fpkm_tracking
    mv genes.fpkm_tracking.sort genes.fpkm_tracking
    mv isoforms.fpkm_tracking.sort isoforms.fpkm_tracking

    # quantification with anaquin Rna Expression
    anaquin RnaExpression -o anaquin_cufflinks -rmix $SEQUINS_ISO_MIX -usequin transcripts.gtf -mix A || (echo "NA" > anaquin_cufflinks/RnaExpression_genes.tsv && echo "NA" > anaquin_cufflinks/RnaExpression_isoforms.tsv && echo "NA" > anaquin_cufflinks/RnaExpression_summary.stats)

echo "FINISH CUFFLINKS: "
date

__CUFF__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# featureCounts
if [ ! -s "feature_counts.txt" ] ; then
    jobid=$(sbatch --export=ALL -J "$fcounts_job" -o "$fcounts_job.o%A" -e "$fcounts_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__FC__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START FEATURECOUNTS: "
date

    FCOUNTS=featureCounts
    FCOUNTS_COMMON="--primary -B -C -p -P --fracOverlap .5 -s 2"
    \$FCOUNTS \$FCOUNTS_COMMON -t 'exon' -g 'gene_id' -a $ANNOTATION -o feature_counts.txt $GENOME_BAM

echo "FINISH FEATURECOUNTS: "
date

__FC__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# kallisto
if [ ! -s "anaquin_kallisto/RnaExpression_summary.stats" ] ; then
    jobid=$(sbatch --export=ALL -J "$kallisto_job" -o "$kallisto_job.o%A" -e "$kallisto_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__KALLISTO__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START KALLISTO: "
date

kallisto quant -i $KALLISTO_INDEX -o kallisto_output $TRIMS_R1 $TRIMS_R2

anaquin RnaExpression -o anaquin_kallisto -rmix $SEQUINS_ISO_MIX -usequin kallisto_output/abundance.tsv -mix A || (echo "NA" > anaquin_kallisto/RnaExpression_genes.tsv && echo "NA" > anaquin_kallisto/RnaExpression_isoforms.tsv && echo "NA" > anaquin_kallisto/RnaExpression_summary.stats)

bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaexp_stats.bash anaquin_kallisto/RnaExpression_summary.stats anaquin_kallisto/RnaExpression_summary.stats.info
bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_neat_comparison.bash anaquin_kallisto/RnaExpression_isoforms.tsv anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv $NEAT_MIX_A

echo "FINISH KALLISTO: "
date

__KALLISTO__
)
    PROCESSING="$PROCESSING,$jobid"
fi

# picard
if [ ! -s "picard.CollectInsertSizes.txt" || ! -s "picard.RnaSeqMetrics.txt" ] ; then
    jobid=$(sbatch --export=ALL -J "$picard_job" -o "$picard_job.o%A" -e "$picard_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__PIC__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START PICARD: "
date

picard CollectInsertSizeMetrics INPUT=Aligned.toGenome.out.bam OUTPUT=picard.CollectInsertSizes.txt HISTOGRAM_FILE=/dev/null
picard CollectRnaSeqMetrics INPUT=Aligned.toGenome.out.bam OUTPUT=picard.RnaSeqMetrics.txt REF_FLAT=$FLAT_REF STRAND="SECOND_READ_TRANSCRIPTION_STRAND"

cat picard.RnaSeqMetrics.txt | grep -A 1 "METRICS CLASS" | sed 1d | tr '\t' '\n' > rna_stats_summary.info
cat picard.RnaSeqMetrics.txt | grep -A 2 "METRICS CLASS" | sed 1d | sed 1d | tr '\t' '\n' | paste rna_stats_summary.info - > tmp.txt && mv tmp.txt rna_stats_summary.info

echo "FINISH PICARD: "
date

__PIC__
)
   PROCESSING="$PROCESSING,$jobid"
fi


# tag counts
if [ ! -s "tagcounts.txt" ] ; then
    jobid=$(sbatch --export=ALL -J "$tagcounts_job" -o "$tagcounts_job.o%A" -e "$tagcounts_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__TC__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START TAGCOUNTS: "
date

if [[ -n $PAIRED ]]; then
    picard MarkDuplicatesWithMateCigar INPUT=Aligned.toGenome.out.bam OUTPUT=/dev/null METRICS_FILE=picard.MarkDuplicatesWithMateCigar.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
else
    picard MarkDuplicates INPUT=Aligned.toGenome.out.bam OUTPUT=/dev/null METRICS_FILE=picard.MarkDuplicatesWithMateCigar.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
fi

python3 $STAMPIPES/scripts/bwa/bamcounts.py Aligned.toGenome.out.bam tagcounts.txt

echo "FINISH TAGCOUNTS: "
date

__TC__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# adapter counts
if [ ! -s "adapter_counts.info" ] ; then
    jobid=$(sbatch --export=ALL -J "$adaptercounts_job" -o "$adaptercounts_job.o%A" -e "$adaptercounts_job.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__

#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: adapter counts"
date

adaptercount=\$(bash "$STAMPIPES/scripts/bam/count_adapters.sh" "Aligned.toGenome.out.bam")
if [ -n \$adapter_count ]; then
        echo -e "adapter\t\$adaptercount" > "adapter_counts.info"
fi

echo "FINISH: adapter counts"
date

__SCRIPT__
)
    PROCESSING="$PROCESSING,$jobid"
fi

# complete
dependencies_full=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterok/g')
sbatch --export=ALL -J "$complete_job" -o "$complete_job.o%A" -e "$complete_job.e%A" $dependencies_full --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__COMPLETE__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START COMPLETE: "
date

    bash $STAMPIPES/scripts/rna-star/aggregate/checkcomplete.sh
#    bash $STAMPIPES/scripts/rna-star/aggregate/upload_counts.bash
    bash $STAMPIPES/scripts/rna-star/aggregate/attachfiles.sh

echo "FINISH COMPLETE: "
date

__COMPLETE__
