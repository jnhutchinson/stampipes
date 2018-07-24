#!/usr/bin/env nextflow

/*
 * This builds all the reference files needed by the pipeline
 */

params.genome = ""
params.read_length = 36
params.outdir = "output"

genome_name = file(params.genome).baseName

process bwa {

  publishDir params.outdir

  input:
  file genome from file(params.genome)

  output:
  file '*'

  script:
  """
  bwa index -p "$genome_name" "$genome"
  """
}

process samtools {
  publishDir params.outdir

  input:
  file genome from file(params.genome)

  output:
  file '*'
  file '*.fai' into fai

  script:
  """
  samtools faidx "$genome"
  """
}

process bowtie_index {

  input:
  file genome from file(params.genome)

  output:
  file "$genome.*" into bowtie_indices

  script:
  """
  bowtie-build "$genome" "$genome"
  """

}

process mappability {

  publishDir params.outdir

  input:
  file genome from file(params.genome)
  val read_length from params.read_length
  file '*' from bowtie_indices

  output:
  file "*bed"
  file "${genome_name}.K${params.read_length}.mappable_only.bed" into mappable

  script:
  """
  awk '/^>/ {\$0=\$1} 1' < "$genome" > cleaned.fa
  perl /hotspot/hotspot-deploy/bin/enumerateUniquelyMappableSpace.pl "$read_length" "$genome" cleaned.fa \
  > "${genome_name}.K${params.read_length}.mappable_only.bed"
  """
}

process density {

  input:
  file genome from file(params.genome)

  output:
  file '*'

  script:
  """
  make all -f "\$STAMPIPES/makefiles/densities/chrombuckets.mk" \
    BWAINDEX="$genome" \
    GENOME="$genome_name" \
    BUCKETS_DIR=.
  """
}

process chrom_sizes {

  input:
  file fai

  output:
  file "${genome_name}.chrom_sizes.bed" into chrom_sizes

  script:
  """
  awk 'BEGIN{OFS="\t"} {print \$1, 0, \$2}' "$fai" > "${genome_name}.chrom_sizes.bed"
  """
}

process hotspot2 {

  container "fwip/hotspot2:latest"

  publishDir params.outdir

  input:
  file chrom_sizes
  file mappable
  val read_length from params.read_length
  val n from 100

  output:
  file '*.starch'

  script:
  """
  extractCenterSites.sh \
    -c "$chrom_sizes" \
    -M "$mappable" \
    -o "${genome_name}.K${read_length}.center_sites.n100.starch" \
    -n "$n"
  """
}
