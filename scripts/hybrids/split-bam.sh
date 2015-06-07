#!/usr/bin/bash
# modules required: samtools, bedops, bedtools, python3
# density_custom.mk must be in the directory this script is run from.
# splits hybrid_bam into two bam files with reads from the genomes
# specified by fa_one and fa_two. Makes bw files for each.

hybrid_bam=$1
fa_one=$2
fa_two=$3

name=`basename $hybrid_bam | rev | cut -d '.' --complement -f1 | rev`
export SAMPLE_NAME=`basename $hybrid_bam | cut -d '.' -f1`
export READLENGTH=36
for ref in $fa_one $fa_two ; do
    genome=`basename $ref | cut -d '.' -f1`
    outbam=$name.$genome.bam
    if ! test -f $ref ; then
        echo "no fasta file found at $ref. Skipping $genome."
    else
        export GENOME=$genome
        export BWAINDEX=$ref
        samtools view $hybrid_bam | awk -v genome=$genome \
            'BEGIN { OFS = "\t" }
            {
                 split($3, ref_chr, "_")
                 ref = ref_chr[1]
                 chr = ref_chr[2]
                 if (ref == genome) {
                     $3 = chr
                     print $0
                 }
             }' | samtools view -S -bt $ref - > $outbam
        make -f density_custom.mk
    fi
done
