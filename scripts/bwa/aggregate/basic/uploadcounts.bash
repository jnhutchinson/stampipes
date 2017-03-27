#!/bin/bash

module load openssl-dev/1.0.1t
module load python/2.7.11

python2 "$STAMPIPES/scripts/lims/upload_aggregation_stats.py" \
  --aggregation "$AGGREGATION_ID" \
  --counts_file "$LIBRARY_NAME.tagcounts.txt"
