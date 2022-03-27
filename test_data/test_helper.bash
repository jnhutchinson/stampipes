#!/bin/bash

function error() {
  echo "$1"
  return 1
}

function require_exe() {
  exitstatus=0
  for exe in "$@" ; do
    if ! command -v "$exe" >/dev/null; then
      exitstatus=1
      echo "Could not find $exe, will not run tests" > /dev/stderr
    fi
  done
  return $exitstatus
}

function cmp_text() {
  name=$1
  echo "Comparing $name..."
  [[ -e "output/$name" ]] || error "output/$name does not exist"
  diff "expected/$name" "output/$name"
}

function cmp_picard() {
  name=$1
  expected=$(grep -v '^#' "expected/$name")
  actual=$(grep -v '^#' "output/$name")

  echo "Comparing $name..."
  diff <(echo "$expected") <(echo "$actual")
}

function cmp_starch() {
  name=$1
  if ! command -v unstarch ; then
    echo "Cannot verify $name, unstarch is not available"
    return 0
  fi
  echo "Comparing $name..."
  cmp <(unstarch "expected/$name") <(unstarch "output/$name") \
    || (echo "$name does not match" ; false)
}

# There are a number of ambiguities that make comparing BAM files naively difficult
function norm_bam() {
  "$samtools" view "$1" |
    awk 'and($2, 4) {$5=0;$6="*"} $2 >= 4096 { $2 -= 4096 } 1' |
    perl -e '$,="\t";' -ane 'print @F[0..10], sort(@F[11..$#F]), "\n"'
}

function cmp_bam() {
  name=$1
  if ! command -v samtools ; then
    echo "Cannot verify $name, samtools is not available"
    return 0
  fi
  echo "Comparing $name..."
  [[ -e "output/$name" ]] || error "output/$name does not exist"
  cmp <(norm_bam "expected/$name") <(norm_bam "output/$name") \
    || (echo "$name does not match" ; false)

}
