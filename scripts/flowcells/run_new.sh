#!/bin/bash

# Where new flowcells are created
flowcell_root=${1:-$SEQUENCER_MOUNT}
# Will only check directories from the last $days_back days
days_back=${2:-10}

set -e -u -o pipefail

if [ -z "$flowcell_root" ] ; then
  echo "Flowcell root not set; set \$SEQUENCER_MOUNT or supply an argument" >&2
  exit 1
fi
if [ ! -d "$flowcell_root" ] ; then
  echo "Flowcell root '$flowcell_root' does not exist" >&2
  exit 1
fi


# Takes an integer and gets the datecode for that many days ago
fmt_day_ago() {
  date --date "$1 days ago" +"%y%m%d"
}

# Takes flowcell_root and days_back
get_recent_flowcells() {
  local flowcell_root=$1
  local days_back=$2
  for day_ago in $(seq "$days_back" -1 0); do
    ls "$flowcell_root/$(fmt_day_ago "$day_ago")_"* -d 2>/dev/null
  done
}

# Takes a directory as argument
flowcell_is_setup() {
  [ -s "$1/run_bcl2fastq.sh" ]
}

# Takes a directory as argument
setup_flowcell() {
  # Run in a subshell to not affect current path
  (
    set -e
    cd "$1"
    bash "$STAMPIPES/scripts/flowcells/setup.sh" -d || exit 1
    nohup bash run_bcl2fastq.sh &
  )
}

# Takes a flowcell directory
log() {
  local fc=$1
  shift;
  echo -e "$(basename "$fc")   \t$*"
}

for fc in $(get_recent_flowcells "$flowcell_root" "$days_back"); do
  if flowcell_is_setup "$fc"; then
    log "$fc" "Already set up"
    continue
  fi

  log "$fc" "Setting up..."
  if setup_flowcell "$fc"; then
    log "$fc" "Setup Successful"
  else
    log "$fc" "Setup Failed"
  fi
done
