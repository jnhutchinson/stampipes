source $MODULELOAD
module load samtools/1.3
module load gcc/4.7.2     # R dependent on this
module load R/3.2.5
module load STAR/2.4.2a   # Just for densities
module load bedops/2.4.19
module load subread/1.5.1 # for featureCounts
module load cufflinks/2.2.1 # for cuffLinks
module load anaquin/2.0

export REFDIR="$(dirname $GENOME_INDEX)"
export STARrefDir="$REFDIR/${STAR_DIR}"
export TARGET_BAM=Aligned.toTranscriptome.out.bam
export GENOME_BAM=Aligned.toGenome.out.bam

sequins_gene_mix="/net/seq/data/genomes/sequins_v1/Sequin_gene_mix.csv"
sequins_isoform_mix="/net/seq/data/genomes/sequins_v1/Sequin_isoform_mix.csv"

if [ -n "$REDO_AGGREGATION" ]; then
    bash $STAMPIPES/scripts/rna-star/aggregate/reset.bash
fi

# create merged fastqs
# placeholder

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

density_job=.AG${AGGREGATION_ID}.star_den
cufflinks_job=.AG${AGGREGATION_ID}.star_cuff
kallisto_job=.AGG${AGGREGATION_ID}.star_kallisto
complete_job=.AGG#${AGGREGATION_ID}.complete
fcounts_job=.AGG${AGGREGATION_ID}.star_fcounts

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

rm -rf "\$TMPDIR"

echo "FINISH DENSITY: "
date

__DEN__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# cufflinks
if [ ! -s "genes.fpkm_tracking" ] ; then
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

    # gene quantification with anaquin Rna Expression
#    ~/anaquin RnaExpression -o anaquin_cufflinks_genes -rmix $sequins_gene_mix -method gene -usequin transcripts.gtf

    # isoform quantification with anaquin Rna Expression
#    ~/anaquin RnaExpression -o anaquin_cufflinks_isoforms -rmix $sequins_isoform_mix -method isoform -usequin transcripts.gtf

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

    # gene quantification with anaquin Rna Expression
#    ~/anaquin RnaExpression -o anaquin_fcounts_genes -rmix $sequins_gene_mix -method gene -usequin #transcripts.gtf 

echo "FINISH FEATURECOUNTS: "
date

__FC__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# kallisto
if [ ! -s "kallisto_output.txt" ] ; then
    jobid=$(sbatch --export=ALL -J "$kallisto_job" -o "$kallisto_job.o%A" -e "$kallisto_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__KALLISTO__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START KALLISTO: "
date

    # run kallisto using gencodev25
    # run kallisto using miitranscriptome

    # gene quantification with anaquin Rna Expression
#    ~/anaquin RnaExpression -o anaquin_fcounts_genes -rmix $sequins_gene_mix -method gene -usequin #transcripts.gtf

echo "FINISH KALLISTO: "
date

__KALLISTO__
)
    PROCESSING="$PROCESSING,$jobid"
fi

dependencies_full=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterok/g')
sbatch --export=ALL -J "$complete_job" -o "$complete_job.o%A" -e "$complete_job.e%A" $dependencies_full --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__COMPLETE__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START COMPLETE: "
date

    bash $STAMPIPES/scripts/rna-star/aggregate/checkcomplete.sh
    bash $STAMPIPES/scripts/rna-star/aggregate/attachfiles.sh

echo "FINISH COMPLETE: "
date

__COMPLETE__
