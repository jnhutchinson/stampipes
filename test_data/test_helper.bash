#!/bin/bash

EXPECTED_DIR=${EXPECTED_DIR:-expected}

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

function get_profile() {
  echo 'test,docker'
}

function compare() {
	local func=$1
	local actual=$2
	local expected=$3

	if [[ ! -f "$actual" ]] ; then echo "# Output file not found: $actual" ; return 1; fi
	if [[ ! -f "$expected" ]] ; then  echo "# Ref file not found: $expected" ; return 1; fi

	cmp <($func "$actual") <($func "$expected") 2>/dev/null
}

function verify() {
  local func=$1
  local filename=$2

  compare "$func" "output/$filename" "$EXPECTED_DIR/$filename"
}

check_starch() {
	unstarch "$1"
}

# There are a number of ambiguities that make comparing BAM files naively difficult
# Some that are addressed here are:
# 1) unsetting undocumented 4096 flag
# 2) Setting mate information for unmapped reads
# 3) Sorting extra tags in consistent order
function check_bam() {
  "$samtools" view "$1" |
    awk 'and($2, 4) {$5=0;$6="*"} $2 >= 4096 { $2 -= 4096 } 1' |
    perl -e '$,="\t";' -ane 'print @F[0..10], sort(@F[11..$#F]), "\n"'
}

check_text() {
	sed '/^#/d' "$1"
}

check_gzipped() {
  zcat < "$1"
}

check_text_sorted() {
	sed '/^#/d;s/\s\+/ /g' "$1" | sort
}
