nextflow.enable.dsl=2

// The number of records in each split file
params.fastq_split_count = 16_000_000

workflow test {
 main:
   // TODO: How to make this work?
   // params.fastq_split_count = 16_000
   meta = [id: "test"]
   r1 = file("test_data/dnase/alignment/r1.fastq.gz")
   r2 = file("test_data/dnase/alignment/r2.fastq.gz")
   input = channel.of([meta, r1, r2])
   split_paired(input).view()
}

def meta_with_read_idx(meta, idx) {
  meta.clone() + [read_idx: idx]
}

workflow split_paired {
  take: // [meta, r1, r2]
    data
  main:
    r1 = data
      .map{ meta, r1, r2 -> [meta_with_read_idx(meta, 1), r1] }
    r2 = data
      .map{ meta, r1, r2 -> [meta_with_read_idx(meta, 2), r2] }
    input = r1.mix(r2)
    input.view { "input_${it}" }

    split_fastq(input)

    output = split_fastq.out
      // Key the output by meta & index, for grouping
      .flatMap { m, files -> {
        i = 1
        r = m.remove("read_idx")
        if (!(files instanceof ArrayList)) {
          files = [files]
        }
        files.collect {[
          [i++, m], // key
            r, // read_idx for sorting
            it, // the file
        ]}
      }}
      // Group together
      .groupTuple(size: 2)
      // Remove extraneous info and emit [meta, r1, r2]
      .map { key, _idx, files -> [ key[1], files.sort { it.name } ] }

  emit: // [meta, r1, r2]
    output
}

process split_fastq {

  tag "${meta?.id}"

  input:
    tuple val(meta), path(fastq)

  output:
    tuple val(meta), file("${name_prefix}*gz")
    
  script:
    fastq_line_count = 4 * params.fastq_split_count
    name_prefix = "${meta?.id}_"
    if (meta.read_idx) {
      name_prefix += "r${meta.read_idx}_"
    }
    """
    zcat "${fastq}" \
    | split -l "${fastq_line_count}" \
      --filter='gzip -1 > \$FILE.gz' \
      - "${name_prefix}"
    """
}
