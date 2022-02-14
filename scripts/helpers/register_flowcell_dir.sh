#!/bin/bash

if [ -n "$1" ] ;then
  cd "$1"
fi

label=$(readlink -f "$PWD" | sed 's/.*FC\([0-9A-Z]*\).*/\1/')
source /home/nelsonjs/code/misc/api_functions.sh

flowcell=$(lims_get_all "flowcell_run/?label=$label" )
id=$(jq .id <<< $flowcell)

python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  --attach_directory "$PWD" \
  --attach_file_objectid "$id" \
  --attach_file_contenttype SequencingData.flowcellrun \
  --attach_file_purpose flowcell-directory
