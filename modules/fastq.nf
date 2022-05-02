nextflow.enable.dsl=2

// The number of records in each split file
params.fastq_split_count = 16_000_000

workflow test {
 main:
   // TODO: How to make this work?
   // params.fastq_split_count = 16_000
   def meta = [id: "test"]
   def r1 = file("test_data/dnase/alignment/r1.fastq.gz")
   def r2 = file("test_data/dnase/alignment/r2.fastq.gz")
   def input = channel.of([meta, r1, r2])
   split_paired(input).view()
}

def meta_with_read_idx(m, idx) {
  def meta_tmp = m.clone()
  meta_tmp.put("read_idx", idx)
  meta_tmp
}

workflow split_paired {
  take: // [meta, r1, r2]
    data
  main:
    def r1 = data
      .map{ meta1, r1fq, r2fq -> [meta_with_read_idx(meta1, 1), r1fq] }
    def r2 = data
      .map{ meta2, r1fq, r2fq -> [meta_with_read_idx(meta2, 2), r2fq] }
    def input = r1.mix(r2)

    //input.view { "Input: $it" }
    split_fastq(input)

    output = split_fastq.out
      // Key the output by meta & index, for grouping
      .flatMap { meta, files ->
        read_idx = meta.read_idx
        m = meta.clone().findAll { it.key != "read_idx" }
        if (!(files instanceof ArrayList)) {
          files = [files]
        }
        files.sort().collect {[
          [it.name.tokenize("_").last(), m], // key
            //read_idx, // read_idx for sorting
            it, // the file
        ]}
      }
      //.view {"Before grouping: ${it[1..-1]}"}
      // Group together
      .groupTuple(size: 2)
      //.view {"After grouping: ${it[1..-1]}"}
      // Remove extraneous info and emit [meta, r1, r2]
      .map { key, files -> [ key[1], files.sort { it.name } ] }

    //output.view { "final output: $it" }

  emit: // [meta, r1, r2]
    output
}

process split_fastq {

  tag "${m?.id}"

  input:
    tuple val(meta), path(fastq)

  output:
    tuple val(m), file("out/*gz")
    
  script:
    // TODO: Why does putting a 'def m' here break this?
    m = meta.clone()
    def fastq_line_count = 4 * params.fastq_split_count
    def name_prefix = "${m?.id}_"
    if (m.read_idx) {
      name_prefix += "r${m.read_idx}_"
    }
    """
    mkdir out
    zcat "${fastq}" \
    | split -l "${fastq_line_count}" \
      --filter='gzip -1 > out/\$FILE.gz' \
      - "${name_prefix}"
    """
}
