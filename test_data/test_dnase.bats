#!/usr/bin/env bats

load test_helper

require_exe nextflow docker

root="$BATS_TEST_DIRNAME/.."
export STAMPIPES=$root

@test 'DNase alignment pipeline' {
  cd $BATS_TEST_DIRNAME/dnase/alignment
  nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile test,docker  -resume

  cmp_picard "MarkDuplicates.picard"
  cmp_picard "CollectInsertSizeMetrics.picard"
  cmp_text "tagcounts.txt"
  cmp_text "subsample.r1.spot.out"
  cmp_starch "density.bed.starch"
  cmp_bam "filtered.bam"
  cmp_bam "marked.bam"

}

@test 'DNase pipeline, single-end' {
  cd $BATS_TEST_DIRNAME/dnase/alignment_singleend
  nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile test,cluster,docker -resume -ansi-log false -process.scratch=false --r2="." --debug

  cmp_picard "MarkDuplicates.picard"
  cmp_picard "tagcounts.txt"
  cmp_picard "subsample.r1.spot.out"
  cmp_starch "density.bed.starch"
  cmp_bam "filtered.bam"
  cmp_bam "marked.bam"

}

@test 'DNase aggregation' {
  cd $BATS_TEST_DIRNAME/dnase/aggregation
  nextflow run "$root/processes/bwa/aggregate/basic.nf" -profile test,docker -resume

  cmp_picard CollectInsertSizeMetrics.picard
  cmp_picard MarkDuplicates.picard

  cmp_starch cutcounts.starch
  cmp_starch density.starch
  cmp_starch fragments.starch

  cmp_bam filtered.bam
  cmp_bam merged.bam

  cmp_text adapter.counts.txt
  cmp_text hs_motifs_svmlight.txt
  cmp_text r1.spot.out
  cmp_text tagcounts.txt

  cmp_text peaks/filtered.SPOT.txt
  cmp_starch peaks/filtered.allcalls.starch
  cmp_starch peaks/filtered.cutcounts.starch
  cmp_starch peaks/filtered.density.starch
  for fdr in 0.05 0.01 0.001 ; do
    cmp_starch peaks/filtered.hotspots.fdr$fdr.starch
    cmp_starch peaks/filtered.peaks.fdr$fdr.narrowpeaks.starch
    cmp_starch peaks/filtered.peaks.fdr$fdr.starch
  done

}
