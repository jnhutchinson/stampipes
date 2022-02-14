#!/usr/bin/env bats

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
export NXF_VER=19.10.0

export STAMPIPES=$root
export HOTSPOT_DIR=/home/solexa/hotspot-hpc/hotspot-distr

@test 'DNase alignment pipeline' {
  cd "$BATS_TEST_DIRNAME/dnase/alignment"
  run nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile "$(get_profile)" -resume
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

  verify check_text "MarkDuplicates.picard"
  verify check_text "CollectInsertSizeMetrics.picard"
  verify check_text "tagcounts.txt"
  verify check_text "subsample.r1.spot.out"
  verify check_starch "density.bed.starch"
  verify check_bam "filtered.cram"
  verify check_bam "marked.cram"

}

@test 'DNase pipeline, single-end' {
  cd $BATS_TEST_DIRNAME/dnase/alignment_singleend
  run nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile "$(get_profile)" -resume -ansi-log false -process.scratch=false --r2="."
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

  verify check_text "MarkDuplicates.picard"
  verify check_text "tagcounts.txt"
  verify check_text "subsample.r1.spot.out"
  verify check_starch "density.bed.starch"
  verify check_bam "filtered.cram"
  verify check_bam "marked.cram"

}

@test 'DNase aggregation' {
  cd "$BATS_TEST_DIRNAME/dnase/aggregation"
  run nextflow run "$root/processes/bwa/aggregate/basic.nf" -profile "$(get_profile)" -resume -ansi-log false
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

  verify check_text CollectInsertSizeMetrics.picard
  verify check_text MarkDuplicates.picard

  verify check_starch cutcounts.starch
  verify check_starch density.starch
  verify check_starch fragments.starch

  verify check_bam filtered.cram
  verify check_bam merged.cram

  verify check_text adapter.counts.txt
  verify check_text hs_motifs_svmlight.txt
  verify check_text r1.spot.out
  verify check_text tagcounts.txt

  verify check_text peaks/nuclear.SPOT.txt
  verify check_starch peaks/nuclear.allcalls.starch
  verify check_starch peaks/nuclear.cutcounts.starch
  verify check_starch peaks/nuclear.density.starch

  verify cmp_text peaks/nuclear.SPOT.txt
  verify cmp_starch peaks/nuclear.allcalls.starch
  verify cmp_starch peaks/nuclear.cutcounts.starch
  verify cmp_starch peaks/nuclear.density.starch
  for fdr in 0.05 0.01 0.001 ; do
    verify check_starch "peaks/nuclear.hotspots.fdr$fdr.starch"
    verify check_starch "peaks/nuclear.peaks.fdr$fdr.starch"
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
