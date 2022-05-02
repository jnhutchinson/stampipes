/// Using the open-source adapter trimmer
/// https://github.com/OpenGene/fastp
process fastp_adapter_trim {

  // TODO: Get a newer version in module system!
  // TODO: Set this module in the RNA pipeline
  //module 'fastp/0.21.0'
  cpus 3

  input:
    tuple val(meta), path(r1), path(r2), val(adapterP5), val(adapterP7)

  output:
    tuple val(meta), path('out.r?.fastq.gz'), emit: fastq
    tuple val(meta), path('fastp.json'),      emit: metrics_json
    tuple val(meta), path('fastp.html'),      emit: metrics_html

  script:
    // TODO: Double-check adapter ordering
    """
    fastp \
      --in1 "${r1}" \
      --in2 "${r2}" \
      --adapter_sequence    "${adapterP7}" \
      --adapter_sequence_r2 "${adapterP5}" \
      --out1 "out.r1.fastq.gz" \
      --out2 "out.r2.fastq.gz" \
      --disable_quality_filtering \
      --disable_length_filtering \
      --thread 3
    """
}

/// Our custom in-house adapter-trimming script
process adapter_trim {
  input:
    tuple val(meta), path(r1), path(r2), val(adapterP5), val(adapterP7)

  output:
    tuple val(meta), path('out.r?.fastq.gz'), emit: fastq
    //path 'out.r2.fastq.gz', emit: trimmed_r2
    tuple val(meta), path('adapter_trimming.txt'), emit: metrics
    tuple val(meta), path('trim.counts.txt'), emit: counts

  script:
  """
  echo -e "P7\t$adapterP7\nP5\t$adapterP5" > adapters.txt

  trim-adapters-illumina \
    -f adapters.txt \
    -1 P5 -2 P7 \
    --threads=3 \
    "$r1" \
    "$r2" \
    out.r1.fastq.gz \
    out.r2.fastq.gz \
    &> adapter_trimming.txt
  awk '{print "adapter-trimmed\t" \$NF * 2}' \
    < adapter_trimming.txt \
    > trim.counts.txt
  """
}

def parse_legacy_adapter_file(adapter_file) {
  // returns two values, p7 and p5
  def adapters = [:]
  adapter_file.readLines().each {
    columns = it.split()
    adapters[columns[0]] = columns[1]
  }
  return [adapters.P7, adapters.P5]
}
