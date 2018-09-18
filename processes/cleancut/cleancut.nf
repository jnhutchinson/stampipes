params.rundir = ""
params.outdir = "output"
params.manifest = "manifest.yaml"
params.specsheet = "spec-sheet.tsv"

scriptdir = "${baseDir}/../../scripts/cleancut"

// Parse spec sheet

Channel
  .fromPath(file(params.specsheet))
  .splitCsv(sep: '\t')
  .subscribe { row ->
    println "${row[0]} ${row[1]} ${row[2]}"
  }

r1_files = Channel.fromPath('i7*R1*.fastq.gz')
r2_files = Channel.fromPath('i7*R2*.fastq.gz')

process merge_pairs {
""" Uses the PEAR utility to merge paired-end reads """
  input:
  file forward from r1_files
  file reverse from r2_files

  output:
  file 'merged.fastq.gz' into merged

  script:
  """
  pear -f "${forward}" -r "${reverse}" -o "merged.fastq.gz"
  """
}


/*
process combine_reads {

  input:
  file script from file("${scriptdir}/preprocess.py")

  output:
  file 'reads.fa' into combined

  script:
  """
  python "${script}" combine-reads \
    "${merged}" \
    "reads.fa" \
    "--unmerged-forward-fq" "${merged_forward}" \
    "--unmerged-reverse-fq" "${merged_reverse}"
  """
}


process align_reads {

  input:
  // TODO
  file amplicon
  file reads from combined

  output:
  // TODO
  file 'alignments.psl'

  script:
  """
  blat -noHead "${amplicon}" "${reads}" alignments.psl
  """
}

process parse_alignments {

  input:
  val mode from 'cinco-de-mayo'
  file reads
  file alignments
  file script from file("${scriptdir}/CC_modes.py")

  output:
  file 'out_aligned.txt' into alignment_results
  file 'out_efficiency.txt' into dimer_efficiency
  file 'out_deletion_profile.txt'
  file 'out_cleavage_profile.txt'

  script:
  """
  python "${script}" \
    "${mode}" \
    "${dimer_id}" \
    "${params.manifest}" \
    "${reads}" \
    "${amplicons}" \
    "${alignments}" \
    --prefix "out"
  """

}


process parse_alignments_cinco_de_mayo {

  when:
  mode == "cinco-de-mayo"

  input:
  file reads
  file alignments

  output:
  file 'out_genotypes.txt'

  script:
  """
  # TODO
  """
}


process gather_efficiencies {

  input:
  file 'efficiency*' from efficiencies.collect()

  output:
  file 'efficiencies.txt'

  script:
  """
  # TODO
  """
}

process gather_results {

  input:
  file 'aligned*' from alignment_results.collect()

  output:
  file 'aligned.txt'

  script:
  """
  # TODO
  """
}
/**/
