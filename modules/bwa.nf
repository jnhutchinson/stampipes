nextflow.enable.dsl=2
params.threads = 2


def get_bwa_index_files(fasta_file) {
  ["", ".amb", ".ann", ".bwt", ".fai", ".pac", ".sa"].collect {
    file("${fasta_file}${it}")
  }
}

workflow BWA_ALIGN {
  // A workflow to automatically pick the "right" BWA alignment strategy
  take:
    metadata
    // Format: meta, r1, r2, genome_fasta_file
    // For single-ended data, r2 should be null or the empty string.
  main:
    metadata
      .map { meta, r1, r2, reference ->
        // Get the references we need
        [meta, r1, r2, reference.name, get_bwa_index_files(reference)]}
      .branch { meta, r1, r2, genome_name, genome_files ->
        // Send different inputs to different BWA methods
        samse: r2 == null || r2 == ""
        sampe: meta.read_length && meta.read_length < 70
        mem: true
    }.set { align_with }

    align_with.sampe | bwa_sampe
    align_with.mem | bwa_mem
    align_with.samse
      .map { meta, r1, r2, genome_name, genome_files ->
        [meta, r1, genome_name, genome_files] }
      | bwa_samse

  emit:
    // Format: [meta, out.bam]
    bwa_samse.out.mix(bwa_sampe.out, bwa_mem.out)
}

process bwa_sampe {

  cpus params.threads

  input:
    tuple val(meta), path(r1), path(r2), val(genome), path('ref/*')

  output:
    tuple val(meta), file('out.bam')

  script:
    """
    bwa aln \
      -Y -l 32 -n 0.04 \
      -t "${params.threads}" \
      "ref/${genome}" \
      "${r1}" \
    > out1.sai

    bwa aln \
      -Y -l 32 -n 0.04 \
      -t "${params.threads}" \
      "ref/${genome}" \
      "${r2}" \
    > out2.sai

    bwa sampe \
      -n 10 -a 750 \
      "ref/${genome}" \
      out1.sai out2.sai \
      "${r1}" "${r2}" \
    | samtools view -t "ref/${genome}.fai" \
      --threads "${params.threads}" \
      -o out.bam
    """
}

process bwa_samse {
  cpus params.threads

  input:
    tuple val(meta), path(r1), val(genome), path('ref/*')

  output:
    tuple val(meta), file('out.bam')

  script:
    """
    bwa aln \
      -Y -l 32 -n 0.04 \
      -t "${params.threads}" \
      "ref/${genome}" \
      "${r1}" \
    > out1.sai

    bwa samse \
      -n 10 \
      "ref/${genome}" \
      out1.sai \
      "${r1}" \
    | samtools view -t "ref/${genome}.fai" \
      --threads "${params.threads}" \
      -o out.bam
    """
}

process bwa_mem {
  cpus params.threads
  input:
    tuple val(meta), path(r1), path(r2), val(genome), path('ref/*')

  output:
    tuple val(meta), file('out.bam')

  script:
    // TODO: -C option for transferring UMI to reads
    // TODO: -R option for setting readgroup
    """
    bwa mem \
      -t "${params.threads}" \
      -M \
      "ref/${genome}" \
      "${r1}" "${r2}" \
    | samtools view -t "ref/${genome}.fai" \
      --threads "${params.threads}" \
      -o out.bam
    """
}
