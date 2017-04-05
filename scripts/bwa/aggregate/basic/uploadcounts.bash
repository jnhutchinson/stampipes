#!/bin/bash
set -e

module load openssl-dev/1.0.1t
module load python/2.7.11

upload() {
  for countsfile in "$@"; do
    python2 "$STAMPIPES/scripts/lims/upload_aggregation_stats.py" \
      --aggregation "$AGGREGATION_ID" \
      --counts_file "$countsfile"
  done
}

# Upload regular stats
upload "$LIBRARY_NAME.tagcounts.txt"

# Upload SPOT2 stats
upload "$LIBRARY_NAME.$GENOME.uniques.sorted.hotspot2.info"

# Upload SPOT1 stats
upload "$LIBRARY_NAME.$GENOME.R1.rand.uniques.sorted.spot.info"

# Upload Adapter stats
upload "$LIBRARY_NAME.adaptercounts.txt"
