nextflow.enable.dsl=2
params.cram_compression_threads = 10
params.cram_version = "3.0"
params.cram_compression_level = 7
params.cram_lossy_names = 0
params.cram_write_index = true
params.cram_compress_other_args = ""  // format:  "key1=value,key2=value2", etc.


workflow test {
  bam = file("test_data/dnase/aggregation/in_1.bam" )
  ref = file("test_data/ref/chr22.fa")
  encode_cram([[], bam, ref])

  encode_cram.out.cram.view { "cram: ${it}" }
  encode_cram.out.cram_with_index.view { "cram with index: ${it}" }
}


process encode_cram {

  module "samtools/1.12"
  cpus Math.ceil(params.cram_compression_threads / 2)

  input:
    tuple val(meta), path(input_bam), path(reference)

  output:
    tuple val(meta), path(output_cram_name), emit: cram
    tuple val(meta), path(output_cram_name), path("${output_cram_name}.crai"), emit: cram_with_index optional true

  script:
    output_cram_name = "${input_bam.baseName}.cram"
    fmt_options = [
      "version=${params.cram_version}",
      "level=${params.cram_compression_level}",
      "lossy_names=${params.cram_lossy_names}",
      "${params.cram_compress_other_args}",
      ].join(",")
    writeindexflag = params.cram_write_index ? "--write-index" : ""
    """
    samtools view \
      -T "${reference}" \
      -C -o "${output_cram_name}" \
      --output-fmt-option "${fmt_options}" \
      --threads "${params.cram_compression_threads}" \
      ${writeindexflag} \
      "${input_bam}"
    """
}
