#!/bin/bash

if [ "$#" -lt "4" ] ; then
    echo "Usage:  $0  input_F.gz input_R.gz output_F.gz output_R.gz" 1>&2
    exit 1
fi
readsF=$1
readsR=$2
outputF=$3
outputR=$4

outputdir=`dirname $outputF`
mkdir -p $outputdir
#if ! [ -d "$outputdir" ] ; then
#    mkdir -p $outputdir
#fi
#
#if [ "$TMPDIR_t" == "" ] ; then
#    TMPDIR_t=/tmp
#fi
#hostname
#date
#echo "TMPDIR_t=$TMPDIR_t"
#
#nameF=`basename $readsF | cut -f1 -d .`
#nameR=`basename $readsR | cut -f1 -d .`
#tmpInputF=$TMPDIR_t/original.$nameF.fq
#tmpInputR=$TMPDIR_t/original.$nameR.fq
#echo "Uncompressing to temp folder" 1>&2
#pwd
#echo -e "zcat -f $readsF > $tmpInputF"
#zcat -f $readsF > $tmpInputF
#ls -altr $readsF $tmpInputF
#zcat -f $readsR > $tmpInputR
#ls -altr $readsR $tmpInputR
#
#if ! [ -s "$tmpInputF" ] ; then
#    echo "Error:  empty $tmpInputF" 1>&2
#    exit 1
#elif ! [ -s "$tmpInputR" ] ; then
#    echo "Error:  empty $tmpInputR" 1>&2
#    exit 1
#fi
#
#tmpOutputF=$TMPDIR_t/trimmed.$nameF.fq
#tmpOutputR=$TMPDIR_t/trimmed.$nameR.fq

adapterFA=$REF_DIR/contamination/IlluminaAdapter_min40bp_revcomp.fa

#$SCRIPT_DIR/fastq-mcf -f -P 33 -p 15 -o $tmpOutputF -o $tmpOutputR $adapterFA $tmpInputF $tmpInputR 
#$SCRIPT_DIR/fastq-mcf -f -P 33 -p 15 -o $outputF -o $outputR $adapterFA <(zcat $readsF) <(zcat $readsR )
$SCRIPT_DIR/fastq-mcf -f -P 33 -p 15 -o $outputF -o $outputR $adapterFA $readsF $readsR

#rm $tmpInputF $tmpInputR 
#if [ -s "$tmpOutputF" ] ; then
#    gzip -9 -c $tmpOutputF > $outputF
#    gzip -9 -c $tmpOutputR > $outputR
#fi
#rm $tmpOutputF $tmpOutputR 
#
#if ! [ -s "$outputF" ] ; then
#    echo "Error:  empty $outputF" 1>&2
#    rm $outputF $outputR
#    exit 1
#fi

#  usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...] 
#  
#  Detects levels of adapter presence, computes likelihoods and
#  locations (start, end) of the adapters.   Removes the adapter
#  sequences from the fastq file, and sends it to stdout.
#  
#  Stats go to stderr, unless -o is specified.
#  
#  If you specify multiple 'paired-end' inputs, then a -o option is
#  required for each.  IE: -o read1.clip.q -o read2.clip.fq
#  
#  Options:
#  	-h      This help
#  	-o FIL  Output file (stats to stdout)
#  	-s N.N  Log scale for clip pct to threshold (2.5)
#  	-t N    % occurance threshold before clipping (0.25)
#  	-m N    Minimum clip length, overrides scaled auto (1)
#  	-p N    Maximum adapter difference percentage (20)
#  	-l N    Minimum remaining sequence length (15)
#  	-L N    Maximum sequence length (none)
#  	-k N    sKew percentage causing trimming (2)
#  	-q N    quality threshold causing trimming (10)
#  	-w N    window-size for quality trimming (1)
#  	-f      force output, even if not much will be done
#  	-0      Set all trimming parameters to zero
#  	-U|u    Force disable/enable illumina PF filtering
#  	-P N    phred-scale (64)
#  	-x N    'N' (Bad read) percentage causing trimming (10)
#  	-R      Don't remove N's from the fronts/ends of reads
#  	-n      Don't clip, just output what would be done
#  	-C N    Number of reads to use for subsampling (100k)
#  	-S FIL  Save clipped reads to file
#  	-d      Output lots of random debugging stuff
#  
#  Increasing the scale makes recognition-lengths longer, a scale
#  of 100 will force full-length recognition.
#  
#  Set the skew (-k) or N-pct (-x) to 100 or 0 to turn it off
