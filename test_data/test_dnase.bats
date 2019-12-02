#!/usr/bin/env bats

load test_helper

root="$BATS_TEST_DIRNAME/.."
export STAMPIPES=$root
export NXF_VER=19.10.0

@test 'DNase alignment pipeline' {
  cd $BATS_TEST_DIRNAME/dnase/alignment
  nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile "$(get_profile)"  -resume -ansi-log false

  verify check_text "MarkDuplicates.picard"
  verify check_text "CollectInsertSizeMetrics.picard"
  verify check_text "tagcounts.txt"
  verify check_text "subsample.r1.spot.out"
  verify check_starch "density.bed.starch"
  verify check_bam "filtered.bam"
  verify check_bam "marked.bam"

}

@test 'DNase pipeline, single-end' {
  cd $BATS_TEST_DIRNAME/dnase/alignment_singleend
  nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile "$(get_profile)" -resume -ansi-log false -process.scratch=false --r2="."

  verify check_text "MarkDuplicates.picard"
  verify check_text "tagcounts.txt"
  verify check_text "subsample.r1.spot.out"
  verify check_starch "density.bed.starch"
  verify check_bam "filtered.bam"
  verify check_bam "marked.bam"

}

@test 'DNase aggregation' {
  cd $BATS_TEST_DIRNAME/dnase/aggregation
  nextflow run "$root/processes/bwa/aggregate/basic.nf" -profile "$(get_profile)"  -resume -ansi-log false

  verify check_text CollectInsertSizeMetrics.picard
  verify check_text MarkDuplicates.picard

  verify check_starch cutcounts.starch
  verify check_starch density.starch
  verify check_starch fragments.starch

  verify check_bam filtered.bam
  verify check_bam merged.bam

  verify check_text adapter.counts.txt
  verify check_text hs_motifs_svmlight.txt
  verify check_text r1.spot.out
  verify check_text tagcounts.txt

  verify check_text peaks/nuclear.SPOT.txt
  verify check_starch peaks/nuclear.allcalls.starch
  verify check_starch peaks/nuclear.cutcounts.starch
  verify check_starch peaks/nuclear.density.starch

  for fdr in 0.05 0.01 0.001 ; do
    verify check_starch peaks/nuclear.hotspots.fdr$fdr.starch
    verify check_starch peaks/nuclear.peaks.fdr$fdr.starch
  done

  # I named these wrong
  verify check_starch peaks/nuclear.peaks.fdr0.001.narrowpeaks.starch
  verify check_starch peaks/nuclear.peaks.fdr0.01.narrowpeaks.starch
  verify check_starch peaks/nuclear.peaks.narrowpeaks.fdr0.05.starch

}

@test 'DNase aggregation single-end' {
  cd $BATS_TEST_DIRNAME/dnase/aggregation_singleend
  nextflow run "$root/processes/bwa/aggregate/basic.nf" -profile "$(get_profile)" --bams "../alignment_singleend/expected/marked.bam" -resume --paired false -process.errorStrategy="ignore" -ansi-log false

  verify check_text MarkDuplicates.picard

  verify check_starch cutcounts.starch
  verify check_starch density.starch
  verify check_starch mm_density.starch
  verify check_starch normalized.density.starch
  verify check_starch normalized.mm_density.starch
  verify check_starch fragments.starch

  verify check_bam filtered.bam
  verify check_bam merged.bam

  verify check_text adapter.counts.txt
  verify check_text hs_motifs_svmlight.txt
  verify check_text r1.spot.out
  verify check_text tagcounts.txt

  verify check_text peaks/nuclear.SPOT.txt
  verify check_starch peaks/nuclear.allcalls.starch
  verify check_starch peaks/nuclear.cutcounts.starch
  verify check_starch peaks/nuclear.density.starch

  for fdr in 0.05 0.01 0.001 ; do
    verify check_starch peaks/nuclear.hotspots.fdr$fdr.starch
    verify check_starch peaks/nuclear.peaks.fdr$fdr.starch
  done

  # I named these wrong
  verify check_starch peaks/nuclear.peaks.fdr0.001.narrowpeaks.starch
  verify check_starch peaks/nuclear.peaks.fdr0.01.narrowpeaks.starch
  verify check_starch peaks/nuclear.peaks.narrowpeaks.fdr0.05.starch
}
