#!/usr/bin/env bats

set -e
here=$BATS_TEST_DIRNAME

@test "Pipeline runs successfully" {
  STAMPIPES=$BATS_TEST_DIRNAME/../../..
  export STAMPIPES
  nextflow run "$STAMPIPES/processes/bwa/aggregate/basic.nf" \
    -profile test,docker \
    -resume
}

# Verify
export EXPECTED_DIR=$BATS_TEST_DIRNAME/expected
load ../../test_helper

@test "Pipeline produces correct output" {
  verify check_text tagcounts.txt
  verify check_text MarkDuplicates.picard
  verify check_text CollectInsertSizeMetrics.picard
  verify check_text r1.spot.out

  verify check_text hs_motifs_svmlight.cols.txt
  verify check_text hs_motifs_svmlight.rows.txt
  verify check_text hs_motifs_svmlight.txt
  verify check_text prox_dist.info

  verify check_bam filtered.bam
  verify check_bam marked.bam

  verify check_starch density.starch
  verify check_starch normalized.density.starch
  verify check_starch mm_density.starch
  verify check_starch normalized.mm_density.starch
  verify check_text preseq.txt
}
