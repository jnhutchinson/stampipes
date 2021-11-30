#!/bin/bash

if [[ $# -ne 1 ]] ; then
  echo "Usage: $0 genome_dir" 1>&2
  exit 1
fi

if ! [[ -d "$1" ]] ; then
  echo "$1 is not a directory" 1>&2
fi

targets=(/home/solexa/stampipes/data $STAMPIPES_DATA)

for target in "${targets[@]}" ; do
  rsync -aLv "$1"/* "$target"
done
