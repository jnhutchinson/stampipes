#!/bin/bash

# TODO: Split this into a makefile? Should be straightforward

if [ "$#" -lt "2" ] ; then
    echo "Usage:  $0 file.bam genome_build"
    exit 1;
fi

bam=$1
# refseq chrom.sizes required by bedGraphToBigWig
refseq=$2

name=`basename $bam | cut -f1 -d .`
destdir=`dirname $bam`
gcovpos=$destdir/coverage.$name.pos.$refseq.bed
gcovneg=$destdir/coverage.$name.neg.$refseq.bed
gcovall=$destdir/coverage.$name.all.$refseq.bed

# First see if this has already been done
#for file in $gcovpos $gcovneg ; do
for file in $gcovall ; do
    subname=${file%%.bed*}
    bigwig=$subname.bigWig
    if [ -s "$bigwig" ] ; then
        echo "Already have $bigwig" 1>&2
        exit
    fi
done

# Create the coverage bed file for each strand
if ! [ -f "$gcovpos" ] ; then
    destpos=$destdir/readspans.$name.pos.$refseq.starch
    destneg=$destdir/readspans.$name.neg.$refseq.starch
    destall=$destdir/readspans.$name.all.$refseq.starch
    if ! [ -f "$destpos" ] ; then
        echo -e "$name $refseq"
        #bamToBed -split -i $bam \
        #    | bbms.pl \
        #    | awk "BEGIN{OFS=\"\\t\"}{ print \$1,\$2,\$3 > \"$destall\"; if (\$6==\"+\") { print \$1,\$2,\$3 > \"$destpos\" } else if (\$6==\"-\") { print \$1,\$2,\$3 > \"$destneg\" }; }"
        $SCRIPT_DIR/splitCoverageByTemplateStrand.pl $bam $destdir/readspans.$name $refseq
        echo "42"
    fi
    echo "44"
    bedops --ec -u $destpos | $SCRIPT_DIR/singleBedFileBaseCoverage.sh | $SCRIPT_DIR/compressBed4.pl > $gcovpos
    bedops --ec -u $destneg | $SCRIPT_DIR/singleBedFileBaseCoverage.sh | $SCRIPT_DIR/compressBed4.pl > $gcovneg
    bedops --ec -u $destall | $SCRIPT_DIR/singleBedFileBaseCoverage.sh | $SCRIPT_DIR/compressBed4.pl > $gcovall
    echo "48"
    if [ -s "$gcovpos" ] ; then
        # cleanup
        rm $destpos $destneg $destall
    fi
fi

chromSizes=$REF_DIR/$refseq/$refseq.chrom.sizes 
if ! [ -f "$chromSizes" ] ; then
    chromSizes=$REF_DIR/$refseq/chrom_sizes.txt 
    if ! [ -f "$chromSizes" ] ; then
        echo "Failed to find $chromSizes " 1>&2
        exit 1
    fi
fi

# Convert to bigWig
for file in $gcovpos $gcovneg $gcovall ; do
    name=${file%%.bed*}
    bigwig=$name.bigWig
    if ! [ -f "$bigwig" ] ; then
        echo "$bigwig"
        bedGraphToBigWig $file $chromSizes $bigwig
    fi
    if [ -s "$bigwig" ] ; then
        # cleanup
        rm $file
    fi
done


