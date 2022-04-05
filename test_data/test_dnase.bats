#!/usr/bin/env bats

export NXF_VER=20.10.0


function get_profile() {
	if [[ "$EXEC_ENV" == "docker" ]] ; then
    echo "test,docker"
  elif [[ "$EXEC_ENV" == "modules" ]] ; then
		echo "test,modules"
  elif command -v docker >/dev/null ; then
		echo >&2 "Using 'docker' for execution environment (override by setting EXEC_ENV)"
    echo "test,docker"
  elif command -v module >/dev/null ; then
		echo >&2 "Using 'modules' for execution environment (override by setting EXEC_ENV)"
    echo "test,modules"
  else
		echo >&2 "Using native execution environment, (override by setting EXEC_ENV)"
    echo "test"
	fi
}

load test_helper

profile=$(get_profile)

require_exe nextflow #docker

root="$BATS_TEST_DIRNAME/.."

export STAMPIPES=$root
export HOTSPOT_DIR=/home/solexa/hotspot-hpc/hotspot-distr

@test 'DNase alignment pipeline' {
  cd "$BATS_TEST_DIRNAME/dnase/alignment"
  rm -rf output/
  run nextflow run "$root/processes/bwa/process_bwa_paired_trimmed.nf" -profile "$profile" -resume -ansi-log false
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

  cmp_picard "MarkDuplicates.picard"
  cmp_picard "CollectInsertSizeMetrics.picard"
  cmp_text "tagcounts.txt"
  cmp_text "subsample.r1.spot.out"
  cmp_starch "density.bed.starch"
  #cmp_bam "filtered.bam"
  #cmp_bam "marked.bam"
  #cmp_bam marked.cram
  cmp_bam filtered.cram

}

@test 'DNase pipeline, single-end' {
  skip "single-end not tested yet"
  cd "$BATS_TEST_DIRNAME/dnase/alignment"
  rm -rf output/
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
  rm -rf output/
  run nextflow run "$root/processes/bwa/aggregate/basic.nf" -profile "$profile" -resume -ansi-log false
  [ "$status" -eq 0 ] || (echo "Output:"; echo "$output" ; false)

  cmp_picard CollectInsertSizeMetrics.picard
  cmp_picard MarkDuplicates.picard

  cmp_starch cutcounts.starch
  cmp_starch density.starch
  cmp_starch fragments.starch

  cmp_bam filtered.cram
  #cmp_bam marked.cram

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
    cmp_starch peaks/nuclear.peaks.fdr$fdr.starch
  done
  cmp_starch peaks/nuclear.peaks.fdr0.001.narrowpeaks.starch
  cmp_starch peaks/nuclear.peaks.fdr0.01.narrowpeaks.starch
  cmp_starch peaks/nuclear.peaks.narrowpeaks.fdr0.05.starch

}
