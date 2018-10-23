#!/usr/bin/env nextflow

params.transcriptbams = ""
params.genomebams = ""
params.outdir = "output"

params.seed = 12345
params.rsem_threads = 8
params.refdir = "/net/seq/data/genomes/human/GRCh38/noalts-sequins/"
starIndexDir = "${params.refdir}/STARgenome-gencode-v25"
rsemIndexDir = "${params.refdir}/RSEMgenome-gencode-v25"
kallistoIndex = "${params.refdir}/ref/kallisto_gencode25_sequins"
annotation = "${params.refdir}/ref/combined.annotation.gtf"
sequinsMix = "/net/seq/data/genomes/sequins_v1/Sequin_isoform_mix.csv"


bams = Channel.create()

num_bams = 0
params.genomebams.tokenize(',').forEach {
  bams << file(it)
  num_bams ++
}
bams.close()

process sort {

  module 'samtools/1.7'
  cpus 8

  input:
  file bam from bams

  output:
  file 'out.bam' into sorted

  script:
  """
  samtools sort \
    "$bam" \
    -o out.bam \
    --threads 8
  """
}


process merge_bams {
  module 'samtools/1.7'

  publishDir params.outdir

  input:
  file 'in*.bam' from sorted.collect()

  output:
  file 'out.bam' into merged

  script:
  if (num_bams == 1) {
    "ln -s in.bam out.bam"
  } else {
    "samtools merge out.bam in*.bam"
  }
}

bam_to_rsem = Channel.create()

to_coverage = Channel.create()
to_kallisto = Channel.create()
to_cufflinks = Channel.create()
to_feature_counts = Channel.create()
to_mark_duplicates = Channel.create()
to_insert_size = Channel.create()
to_adapter_counts = Channel.create()
merged.into(to_coverage, to_kallisto, to_cufflinks, to_feature_counts, to_mark_duplicates, to_insert_size, to_adapter_counts)

process coverage {
  module 'STAR', 'samtools/1.7'

  publishDir params.outdir

  input:
  file bam from to_coverage

  output:
  file 'Signal*bg' into bedgraph

  script:
    """
    samtools index "$bam"
    STAR \
      --runMode inputAlignmentsFromBAM \
      --inputBAMfile "$bam" \
      --outWigType bedGraph --outWigStrand Stranded \
      --outFileNamePrefix ./ \
      --outWigReferencesPrefix chr

#--outTmpDir STAR
    """
}

to_starch = Channel.create()
to_bigwig = Channel.create()

bedgraph.flatten().into(to_starch, to_bigwig)

process coverage_starch {

  module 'bedops/2.4.19'
  publishDir params.outdir

  input:
  file bg from to_starch

  output:
  file '*starch'

  script:
  """
  sort-bed - < "$bg" \
  | starch - > "${bg}.starch"
  """
}

process coverage_density {

  module 'kentutil/302'
  publishDir params.outdir

  input:
  file bg from to_bigwig
  file cnl from file("${starIndexDir}/chrNameLength.txt")
  output:
  file '*bw'

  script:
  """
  grep '^chr' "$bg" | sort -k1,1 -k2,2n > chrom.bg
  grep '^chr' "$cnl" > chrNL.txt
  bedGraphToBigWig "chrom.bg" chrNL.txt "${bg}.bw"
  """
}

process trimmed_fastq {

  module 'samtools/1.7'

  input:
  file bam from to_kallisto

  output:
  set file('r1.fastq'), file('r2.fastq') into trimmed_fastq

  script:
  """
  samtools sort -n "$bam" --threads 8 \
  | tee >(
    samtools view -u -f  64 "$bam" | samtools fastq - > r1.fastq
      ) \
  | samtools view -u -f 128 "$bam" | samtools fastq - > r2.fastq
  """
}

for_kallisto = Channel.create()
for_ribosomal = Channel.create()
trimmed_fastq.into(for_kallisto, for_ribosomal)

process kallisto {

  module 'kallisto/0.43.1'
  publishDir params.outdir

  input:
  set file(r1), file(r2) from for_kallisto
  file index from file(kallistoIndex)

  output:
  file 'kallisto_output/*'

  script:
  """
  kallisto quant -i \
    "$index" \
    -o kallisto_output \
    "$r1" "$r2"
  """
}

process insert_size {

  module 'jdk/1.8.0_92', 'picard/2.8.1', 'R/3.2.5'
  publishDir params.outdir

  input:
  file bam from to_insert_size

  output:
  file 'picard.CollectInsertSizes.txt'

  script:
  """
  picard CollectInsertSizeMetrics \
  INPUT="$bam" \
  OUTPUT=picard.CollectInsertSizes.txt \
  HISTOGRAM_FILE=/dev/null
  """
}

process adapter_counts {

  publishDir params.outdir
  input:
  file bam from to_adapter_counts

  output:
  file 'adapter_counts.info'

  script:
  """
  export STAMPIPES=$baseDir/../../..
  echo -ne "adapter\t" > adapter_counts.info
  bash "$baseDir/../../../scripts/bam/count_adapters.sh" "$bam" >> adapter_counts.info
  """

}

process cufflinks {

  module 'cufflinks/2.2.1'
  publishDir params.outdir

  input:
  file input from to_cufflinks
  file annotation from file(annotation)

  output:
  file 'transcripts.gtf' into transcripts
  file 'skipped.gtf'
  file 'genes.fpkm_tracking'
  file 'isoforms.fpkm_tracking'


  script:
  """
  cufflinks \
    --no-update-check \
    --library-type fr-firststrand \
    --GTF "${annotation}" \
    --num-threads 8 \
    "${input}"
  """
}

process anaquin {

  module 'anaquin/2.0.1'
  publishDir params.outdir

  input:
  file transcripts
  file sequins_mix from file(sequinsMix)

  output:
  file 'anaquin_cufflinks/RnaExpression*'

  script:
  """
  anaquin RnaExpression \
    -o anaquin_cufflinks \
    -rmix "${sequins_mix}" \
    -usequin "${transcripts}" \
    -mix A
  """
}

process ribosomal {
  module 'bowtie/1.0.0'
  publishDir params.outdir

  input:
  file 'index.1.ebwt' from file('/net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA.1.ebwt')
  file 'index.2.ebwt' from file('/net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA.2.ebwt')
  file 'index.3.ebwt' from file('/net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA.3.ebwt')
  file 'index.4.ebwt' from file('/net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA.4.ebwt')
  file 'index.rev.1.ebwt' from file('/net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA.rev.1.ebwt')
  file 'index.rev.2.ebwt' from file('/net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA.rev.2.ebwt')
  set file(r1), file(r2) from for_ribosomal
  val seed from 1234567

  output:
  file 'ribosomal_counts.info'

  script:
  """
  sleep 20;
  echo -ne "ribosomal-RNA\t" > ribosomal_counts.info
  bowtie \
    -n 3 -e 140 \
    --threads 1 \
    --seed "$seed" \
    index \
    -1 "$r1" \
    -2 "$r2" \
    | wc -l \
    >> ribosomal_counts.info
  """

}

process mark_duplicates {
  module 'jdk/1.8.0_92', 'picard/2.8.1'
  publishDir params.outdir

  input:
  file bam from to_mark_duplicates

  output:
  file 'picard.MarkDuplicates.txt'

  script:
  """
  picard MarkDuplicates \
    INPUT="$bam" \
    OUTPUT=/dev/null \
    METRICS_FILE=picard.MarkDuplicates.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
  """

}

process feature_counts {

  module 'subread/1.5.1'
  publishDir params.outdir

  input:
  file input from to_feature_counts
  file annotation from file(annotation)

  output:
  file 'feature_counts.txt'

  script:
  """
  featureCounts \
    --primary \
    -B \
    -C \
    -p \
    -P \
    --fracOverlap .5 \
    -s 2 \
    -t exon \
    -g gene_id \
    -a "${annotation}" \
    -o feature_counts.txt \
    "${input}"
  """
}
