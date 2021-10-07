#!/bin/bash

image=quay.io/biocontainers/adapterremoval:2.3.1--hf3e9acd_0

if [[ $# -lt 2 ]] ; then
  echo "Usage: $0 r1.fastq.gz r2.fastq.gz [adapter1-seq adapter2-seq]" 
  exit 1
fi

if ! command -v docker &>/dev/null ; then
  echo "You need docker; try ssh'ing to another node, like hpcA10" 
  exit 1
fi

f1=$(readlink -f "$1")
f2=$(readlink -f "$2")
a1=$3
a2=$4

f1dir=$(dirname "$f1")
f2dir=$(dirname "$f2")
aparams=()
if [[ -n "$a1" && -n "$a2" ]] ; then
  aparams=(--adapter1 $a1 --adapter2 $a2)
fi

docker run -it --rm \
  -v "$f1dir:$f1dir" \
  -v "$f2dir:$f2dir" \
  "$image" \
  AdapterRemoval --identify-adapters --file1 "$f1" --file2 "$f2" \
  "${aparams[@]}"

