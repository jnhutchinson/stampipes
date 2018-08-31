#!/usr/bin/env nextflow

/*
 * This builds all the reference files needed by the pipeline
 */

params.genome = ""
params.readlength = 36
params.outdir = "output"

genome_name = file(params.genome).baseName
readlengths = params.readlength.tokenize(',')

process bwa {

  publishDir "${params.outdir}/bwa"

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
  publishDir "${params.outdir}/samtools"

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

  publishDir "${params.outdir}/annotations"

  input:
  file genome from file(params.genome)
  val read_length from Channel.from(readlengths)
  file '*' from bowtie_indices

  output:
  file "*bed"
  set val(read_length), file("${genome_name}.K${read_length}.mappable_only.bed") into mappable

  script:
  """
  awk '/^>/ {\$0=\$1} 1' < "$genome" > cleaned.fa
  perl "\$HOTSPOT_DIR/hotspot-deploy/bin/enumerateUniquelyMappableSpace.pl" "$read_length" "$genome" cleaned.fa \
    | awk -f \$STAMPIPES/awk/merge_adjacent_bed.awk \
    | sed 's/\\s\$//' \
    | sort-bed - \
    > "${genome_name}.K${read_length}.mappable_only.bed"
  """
}

process density {

  publishDir "${params.outdir}/densities"

  input:
  file fai

  output:
  file '*'

  script:
  """
  make all -f "\$STAMPIPES/makefiles/densities/chrombuckets.mk" \
    FAI="$fai" \
    GENOME="$genome_name" \
    BUCKETS_DIR=.
  """
}

process chrom_sizes {

  publishDir "${params.outdir}/hotspot2"

  input:
  file fai

  output:
  file "${genome_name}.chrom_sizes.bed" into chrom_sizes

  script:
  """
  awk 'BEGIN{OFS="\t"} {print \$1, 0, \$2}' "$fai" | sort-bed - > "${genome_name}.chrom_sizes.bed"
  """
}

process chromInfo {

  publishDir "${params.outdir}/annotations"

  input:
  file chrom_sizes

  output:
  file "${genome_name}.chromInfo.bed"

  script:
  """
  awk 'BEGIN{OFS="\t"} {print \$1, \$2, \$3, \$1}' "$chrom_sizes" > "${genome_name}.chromInfo.bed"
  """
}

process hotspot2 {

  container "fwip/hotspot2:latest"

  publishDir "${params.outdir}/hotspot2"

  input:
  file chrom_sizes
  set val(read_length), file(mappable) from mappable
  val n from 100

  output:
  file '*.starch'

  script:
  """
  extractCenterSites.sh \
    -c "$chrom_sizes" \
    -M "$mappable" \
    -o "${genome_name}.K${read_length}.center_sites.n${n}.starch" \
    -n "$n"
  """
}
