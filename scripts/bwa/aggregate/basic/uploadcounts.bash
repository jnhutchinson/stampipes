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
if [[ -n "$PAIRED" ]]; then
        upload "$LIBRARY_NAME.$GENOME.R1.rand.uniques.sorted.spot.info"
else
        upload "$LIBRARY_NAME.$GENOME.rand.uniques.sorted.spot.info"
fi

# Upload Adapter stats
upload "$LIBRARY_NAME.adaptercounts.txt"

# Upload Preseq stats
if [ -e "$LIBRARY_NAME.uniques.preseq.targets.txt" ]; then
    upload "$LIBRARY_NAME.uniques.preseq.targets.txt"
fi

# Upload Insert Size stats
if [ -e $LIBRARY_NAME.CollectInsertSizeMetrics.picard.info ]; then
    upload "$LIBRARY_NAME.CollectInsertSizeMetrics.picard.info"
fi

# Upload Proximal / Distal
if [ -e $LIBRARY_NAME.$GENOME.uniques.sorted.proxdist.txt ]; then
    upload "$LIBRARY_NAME.$GENOME.uniques.sorted.proxdist.txt"
fi
