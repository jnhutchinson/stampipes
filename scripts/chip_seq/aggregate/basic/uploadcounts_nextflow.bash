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

maybe_upload() {
  for countsfile in "$@"; do
    if [ -e "$countsfile" ] ; then
      upload "$countsfile"
    fi
  done
}

# Upload regular stats
upload "tagcounts.txt"

# Upload macs stats
# TODO

# Upload SPOT1 stats
upload r1.hotspot.info

# Upload Adapter stats
upload "adapter.counts.txt"

# Upload Preseq stats
maybe_upload "preseq_targets.txt"

# Upload Insert Size stats
maybe_upload "CollectInsertSizeMetrics.picard.info"

# Upload Proximal / Distal
maybe_upload prox_dist.info
