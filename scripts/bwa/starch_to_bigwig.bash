#!/bin/bash

set -e -u -o pipefail

starch=$1
bigwig=$2
fai=$3
bin=${4:-}

set -x
if [[ -z "$bin" ]] ; then
  unstarch "$starch" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5}' \
    | grep -v chrM \
    > tmp.wig
else
  unstarch "$starch" \
    | awk -v "binI=$bin" -f "$STAMPIPES/awk/bedToWig.awk" \
    > tmp.wig
fi

wigToBigWig -clip tmp.wig "${fai}" "${bigwig}"

rm tmp.wig
