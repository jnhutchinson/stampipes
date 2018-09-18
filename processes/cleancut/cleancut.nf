params.rundir = ""
params.outdir = "output"
params.manifest = "manifest.yaml"
params.specsheet = "spec-sheet.tsv"
params.minreads = 1000

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
  set file('merged.assembled.fastq'), file('merged.unassembled.forward.fastq'), file('merged.unassembled.reverse.fastq') into merged

  script:
  """
  pear -f "${forward}" -r "${reverse}" -o "merged"
  """
}


process combine_reads {

  input:
  file script from file("${scriptdir}/preprocess.py")
  set file(merged), file(unassembled_forward), file(unassembled_reverse) from merged

  output:
  file 'reads.fa' into combined

  script:
  """
  python "${script}" combine-reads \
    "${merged}" \
    "reads.fa" \
    "--unmerged-forward-fq" "${unassembled_forward}" \
    "--unmerged-reverse-fq" "${unassembled_reverse}"
  """
}

process align_reads {

  input:
  file amplicon from file('amplicons.fa')
  file reads from combined

  output:
  // TODO
  set file(reads), file('alignments.psl') into to_parse

  script:
  """
  blat -noHead "${amplicon}" "${reads}" alignments.psl
  """
}

methods = ['dimer-efficiency', 'cinco-de-mayo']
process parse_alignments {

  input:
  //val mode from 'cinco-de-mayo'
  each mode from methods
  val dimer_id from 'fake'
  file amplicon from file('amplicons.fa')
  file script from file("${scriptdir}/CC_modes.py")
  file manifest from file(params.manifest)

  set file(reads), file(alignments) from to_parse

  output:
  file 'out_aligned.txt' optional true into alignment_results
  file 'out_efficiency.txt' optional true into dimer_efficiencies
  file 'out_deletion_profile.txt' optional true
  file 'out_cleavage_profile.txt' optional true
  file 'out_genotypes.txt' optional true into genotype_results

  script:
  """
  python "${script}" \
    "${mode}" \
    "${dimer_id}" \
    "${manifest}" \
    "${reads}" \
    "${amplicon}" \
    "${alignments}" \
    --prefix "out"
  """

}

/*
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
*/


process gather_efficiencies {

  publishDir params.outdir

  input:
  val minreads from params.minreads
  file 'efficiency*' from dimer_efficiencies.collect()

  output:
  file 'efficiencies.txt'

  script:
  """
  cat efficiency* \
  | sed '/^dimer/d' \
  | awk 'BEGIN{
    OFS="\t"
    print "i7 index\ti5 index\tdimer(s)\tindel\tHDR\tmicro\ttotal\tindel %\tHDR %\tmicro %"
  }
    \$5 > $minreads {
      print "ID1", "ID2", \$0
    }' \
  > efficiencies.txt

  """
}

process gather_results {

  publishDir params.outdir

  input:
  file 'aligned*' from alignment_results.collect()

  output:
  file 'aligned.txt'

  script:
  """
  cat aligned* \
  | awk 'BEGIN{OFS="\t"} {print "i7", "i5", "dimer", \$0}' \
  > aligned.txt
  """
}
/**/

process gather_genotypes {
  publishDir params.outdir

  input:
  file 'genotypes*' from genotype_results.collect()

  output:
  file 'genotypes.txt'

  script:
  """
  cat genotypes* \
  | sed '/^dimer/d' \
  | awk 'BEGIN {OFS="\t"} {print "i7", "i5", "dimer", \$0}' \
  > genotypes.txt
  """
}
