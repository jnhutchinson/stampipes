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

# upload stats
upload "tagcounts.txt"
upload "adapter_counts.info"
upload "rna_stats_summary.info"

# upload sequins stats
if [[ -n "$SEQUINS_REF" ]]; then
    upload "anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info"
    upload "anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info"
    upload "anaquin_star/RnaAlign_summary.stats.info"
fi