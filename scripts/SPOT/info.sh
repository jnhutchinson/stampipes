#!/bin/bash

starch=$1
prefix=$2
SPOT=$3
if [[ -z "$starch" || "$(unstarch --is-starch "$starch")" == "0" ]]; then
  echo "Couldn't read starch file '$starch'" >&2
  exit 1
fi

echo -e "$prefix-num-bases\t$(unstarch --bases    "$starch")"
echo -e "$prefix-num-spots\t$(unstarch --elements "$starch")"

if [[ -n "$SPOT" ]]; then
  echo -e "$prefix-SPOT\t$SPOT"
fi
