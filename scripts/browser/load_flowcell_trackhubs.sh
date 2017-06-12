#!/bin/bash

set -e -o pipefail

flowcell_name=$1
priority='400' # the concept of priority is somewhat depreciated with trackhubs
config=$STAMPIPES/config/ucsc_browser/trackhub_flowcells.config	# only one config file for the foreseeable future

flowcell_dir=$(ls $FLOWCELLS/FC${flowcell_name}_*tag -d)

die(){
  echo $@ >&2
  exit 1
}

if [[ -z "$priority" ]] ; then
  die "Usage: $0 flowcell_label priority"
fi

if [[ -z "$flowcell_dir" ]] ; then
  die "Couldn't find flowcell dir for $flowcell"
fi

source "$PYTHON3_ACTIVATE"

cd "$flowcell_dir"

python "$STAMPIPES/scripts/lims/get_processing.py" -f "$flowcell_name" --quiet
python "$STAMPIPES/scripts/browser/make_trackhubs_for_flowcell.py" -p "$priority" -c "$config" --quiet
bash "$STAMPIPES/scripts/browser/make_flowcell_hub.sh" "$flowcell_name"