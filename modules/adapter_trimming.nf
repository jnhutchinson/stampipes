/// Our custom in-house adapter-trimming script
process adapter_trim {
  input:
    tuple path(r1), path(r2), val(adapterP5), val(adapterP7)

  output:
    path 'out.r?.fastq.gz', emit: fastq
    //path 'out.r2.fastq.gz', emit: trimmed_r2
    path 'adapter_trimming.txt', emit: metrics

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
  """
}
