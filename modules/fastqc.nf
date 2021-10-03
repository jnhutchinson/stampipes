#!/usr/bin/env nextflow


params.fastq = "in.fq.gz"
params.fastqc = "out.fastqc.zip"
params.umi_report = "out.umi.gz"

fastq_file = file(params.fastq)
fastqc_output = "${fastq_file.baseName.take(fastq_file.baseName.indexOf('.'))}_fastqc.zip"

STAMPIPES = "/home/nelsonjs/code/stampipes-nextflow"

process fastqc {

  publishDir '.', mode: 'move'
  scratch true

  container = '/home/nelsonjs/tmp/sing/fastqc.img'

//  module 'jdk/1.8.0_92'
//  module 'picard/2.8.1'
//  module 'fastqc/0.11.5'

  input:
  file fastq_file

  output:
  file "${params.fastqc}"

  """
  time fastqc -t 8 --noextract --nogroup --casava "$fastq_file"
  mv "${fastqc_output}" "${params.fastqc}"
  """
}

process umi {
  publishDir '.', mode: 'move'

  input:
  file fastq_file

  output:
  file "${params.umi_report}"


  """
  zcat "${params.fastq}" \
    | grep "^@" \
    | cut -f 2 -d "+" \
    | sort \
    | uniq -c \
    | sort -n -r \
    | gzip -c \
    > "${params.umi_report}"
  """
}
