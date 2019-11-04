#!/bin/bash

function compare() {
	local func=$1
	local actual=$2
	local expected=$3

	if [[ ! -f "$actual" ]] ; then echo "Nonexistent: $actual" ; return 1; fi
	if [[ ! -f "$expected" ]] ; then  echo "Nonexistent: $expected" ; return 1; fi

	cmp <($func "$actual") <($func "$expected") 2>/dev/null || echo "$func: $2 $3 differ"
}

function verify() {
  local func=$1
  local filename=$2

  compare "$func" "output/$filename" "expected/$filename"
}

check_starch() {
	unstarch "$1"
}

check_bam() {
	samtools view "$1"
}

check_text() {
	sed '/^#/d' "$1"
}

check_text_sorted() {
	sed '/^#/d;s/\s\+/ /g' "$1" | sort
}
