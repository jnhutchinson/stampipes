#!/bin/bash

EXPECTED_DIR=${EXPECTED_DIR:-expected}

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

check_bam() {
	samtools view "$1" \
        | awk '$2 >= 4096 { $2 -= 4096 } 1' \
        |  perl -ne '@x=split/\s/; print join " ",sort @x, "\n"'
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
