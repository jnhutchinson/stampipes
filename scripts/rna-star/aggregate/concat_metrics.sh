# makes a concactenated version of all relevant rna seq metrics
# shellcheck disable=SC2002  # It's okay to have a useless use of cat

rm -f metrics.info
{
    echo -ne "fastq-total-reads\t"
    zcat -f r1.fq.gz | wc -l | awk '{print 2*$1/4}'

    echo -ne "picard2-percent-duplication\t"
    cat picard.MarkDuplicates.txt | grep 'READ_PAIR' -A 1 | sed 1d | cut -f 9

    echo -ne "fcounts-assigned\t"
    cat feature_counts.txt.summary | grep 'Assigned' | cut -f 2 | awk '{print $1}'

    echo -ne "kallisto-palign-default\t"
    cat kallisto.log | grep 'reads pseudoaligned' | tr ' ' '\t' | cut -f 5 | sed -e 's/,//g' | awk '{print $1*2}'

    cat ribosomal_counts.info | grep 'ribosomal'
    cat adapter_counts.info | grep 'adapter'

    if [ -s anaquin_star/RnaAlign_summary.stats.info ]; then
        cat anaquin_star/RnaAlign_summary.stats.info | grep 'sequins-dilution'
        cat anaquin_star/RnaAlign_summary.stats.info | grep 'sequins-base-level-sensitivity'
        cat anaquin_star/RnaAlign_summary.stats.info | grep 'sequins-base-level-precision'
    fi

    if [ -s anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info ]; then
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv.info | grep 'neat-mixA-mean-spearman'
    fi

    if [ -s anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info ]; then
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-isoforms-log2-pearson-cor'
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-genes-slope'
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-genes-log2-pearson-cor'
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-percent-isoforms-found'
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-percent-genes-found'
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-detection-sensitivity-isoforms'
        cat anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info | grep 'sequins-detection-sensitivity-genes'
    fi
} > metrics.info
