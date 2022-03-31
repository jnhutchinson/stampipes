nextflow.enable.dsl=2
params.threads = 2


def get_bwa_index_files(fasta_file) {
  ["", ".amb", ".ann", ".bwt", ".fai", ".pac", ".sa"].collect {
    file("${fasta_file}${it}")
  }
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
    | samtools view -b -t "ref/${genome}.fai" - \
    > out.bam
    """
}
