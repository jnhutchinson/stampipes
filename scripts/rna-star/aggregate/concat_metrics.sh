# makes a concactenated version of all relevant rna seq metrics

rm -f metrics.info
echo -ne "fastq-total-reads\t" > metrics.info
zcat -f r1.fq.gz | wc -l | awk '{print 2*$1/4}' >> metrics.info

echo -ne "picard2-percent-duplication\t" >> metrics.info
cat picard.MarkDuplicates.txt | grep 'READ_PAIR' -A 1 | sed 1d | cut -f 9 >> metrics.info

echo -ne "fcounts-assigned\t" >> metrics.info
cat feature_counts.txt.summary | grep 'Assigned' | cut -f 2 | awk '{print $1}' >> metrics.info

echo -ne "kallisto-palign-default\t" >> metrics.info
cat kallisto.log | grep 'reads pseudoaligned' | tr ' ' '\t' | cut -f 5 | sed -e 's/,//g' | awk '{print $1*2}' >> metrics.info

cat ribosomal_counts.info | grep 'ribosomal' >> metrics.info
cat adapter_counts.info | grep 'adapter' >> metrics.info

if [ -s anaquin_star/RnaAlign_summary.stats.info ]; then
cat anaquin_star/RnaAlign_summary.stats.info | grep 'sequins-dilution' >> metrics.info
cat anaquin_star/RnaAlign_summary.stats.info | grep 'sequins-base-level-sensitivity' >> metrics.info
cat anaquin_star/RnaAlign_summary.stats.info | grep 'sequins-base-level-precision' >> metrics.info
fi

if [ -s anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info ]; then
cat anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info | grep 'neat-mixA-mean-spearman' >> metrics.info
fi

if [ -s anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info ]; then
cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-isoforms-log2-pearson-cor' >> metrics.info
cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-genes-slope' >> metrics.info
cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-genes-log2-pearson-cor' >> metrics.info
cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-percent-isoforms-found' >> metrics.info
cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-percent-genes-found' >> metrics.info
cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-detection-sensitivity-isoforms' >> metrics.info
cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-detection-sensitivity-genes' >> metrics.info
fi
