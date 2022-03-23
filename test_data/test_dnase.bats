#!/usr/bin/env bats

export NXF_VER=20.10.0

load test_helper

if command -v docker >/dev/null ; then
  profile=test,docker
elif command -v module >/dev/null ; then
  profile=test,modules
  module purge
  module load jdk nextflow samtools/1.12 bedops
else
  echo "Neither docker nor modules available, can't test."
  exit 1
fi

require_exe nextflow #docker

root="$BATS_TEST_DIRNAME/.."

export STAMPIPES=$root
export HOTSPOT_DIR=/home/solexa/hotspot-hpc/hotspot-distr

@test 'DNase alignment pipeline' {
  cd "$BATS_TEST_DIRNAME/dnase/alignment"
  run nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile "$profile" -resume -ansi-log false
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

  cmp_picard "MarkDuplicates.picard"
  cmp_picard "CollectInsertSizeMetrics.picard"
  cmp_text "tagcounts.txt"
  cmp_text "subsample.r1.spot.out"
  cmp_starch "density.bed.starch"
  #cmp_bam "filtered.bam"
  #cmp_bam "marked.bam"
  cmp_bam marked.cram

}

@test 'DNase pipeline, single-end' {
  skip "single-end not tested yet"
  cd "$BATS_TEST_DIRNAME/dnase/alignment"
  run nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile "$profile" --r2 "" -resume -ansi-log false
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

#  cmp_picard "MarkDuplicates.picard"
#  cmp_picard "CollectInsertSizeMetrics.picard"
#  cmp_picard "tagcounts.txt"
#  cmp_picard "subsample.spot.out"
#  cmp_starch "density.bed.starch"
#  cmp_bam "filtered.bam"
#  cmp_bam "marked.bam"

}

@test 'DNase aggregation' {
  cd "$BATS_TEST_DIRNAME/dnase/aggregation"
  run nextflow run "$root/processes/bwa/aggregate/basic.nf" -profile "$profile" -resume -ansi-log false
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

  cmp_picard CollectInsertSizeMetrics.picard
  cmp_picard MarkDuplicates.picard

  cmp_starch cutcounts.starch
  cmp_starch density.starch
  cmp_starch fragments.starch

  cmp_bam filtered.cram
  cmp_bam merged.cram

  cmp_text adapter.counts.txt
  cmp_text hs_motifs_svmlight.txt
  cmp_text r1.spot.out
  cmp_text tagcounts.txt

  cmp_text peaks/nuclear.SPOT.txt
  cmp_starch peaks/nuclear.allcalls.starch
  cmp_starch peaks/nuclear.cutcounts.starch
  cmp_starch peaks/nuclear.density.starch
  for fdr in 0.05 0.01 0.001 ; do
    cmp_starch peaks/nuclear.hotspots.fdr$fdr.starch
    cmp_starch peaks/nuclear.peaks.fdr$fdr.narrowpeaks.starch
    cmp_starch peaks/nuclear.peaks.fdr$fdr.starch
  done

}
