#!/usr/bin/env bats

load test_helper

require_exe nextflow docker

root="$BATS_TEST_DIRNAME/.."

@test 'DNase alignment pipeline' {
  cd $BATS_TEST_DIRNAME/dnase/alignment
  run nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile test,docker -resume
  [ "$status" -eq 0 ]

  cmp_picard "MarkDuplicates.picard"
  cmp_picard "CollectInsertSizeMetrics.picard"
  cmp_text "tagcounts.txt"
  cmp_text "subsample.r1.spot.out"
  cmp_starch "density.bed.starch"
  cmp_bam "filtered.bam"
  cmp_bam "marked.bam"

}

@test 'DNase pipeline, single-end' {
  skip "single-end not tested yet"
  cd $BATS_TEST_DIRNAME/dnase/alignment
  run nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile test,docker --r2 "" -resume
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
  cd $BATS_TEST_DIRNAME/dnase/aggregation
  run nextflow run "$root/processes/bwa/aggregate/basic.nf" -profile test,docker -resume

  cmp_picard CollectInsertSizeMetrics.picard
  cmp_picard MarkDuplicates.picard

  cmp_starch cutcounts.starch
  cmp_starch density.starch
  cmp_starch fragments.starch

  cmp_bam filtered.bam
  cmp_bam merged.bam

  cmp_text adapter.counts.txt
  cmp_text hs_motifs_svmlight.txt
  cmp_text subsample.spot.out
  cmp_text tagcounts.txt

  cmp_text peaks/filtered.SPOT.txt
  cmp_starch peaks/filtered.allcalls.starch
  cmp_starch peaks/filtered.cutcounts.starch
  cmp_starch peaks/filtered.density.starch
  cmp_starch peaks/filtered.hotspots.fdr0.05.starch
  cmp_starch peaks/filtered.peaks.narrowpeaks.starch
  cmp_starch peaks/filtered.peaks.starch

}
