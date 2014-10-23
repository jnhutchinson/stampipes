#FINAL_BAM=${SAMPLE_NAME}.tophat.bam

export SCRIPT_DIR="$STAMPIPES/scripts/tophat"
export REF_DIR="$STAMPIPES/data/tophat/refseq"

if [ ! -e ${SAMPLE_NAME}_R1_fastqc -o ! -e ${SAMPLE_NAME}_R2_fastqc ]; then
qsub -N ".fq${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/fastqc.mk
__SCRIPT__
fi

nohup qmake -cwd -V -q all.q -N .th-$SAMPLE_NAME -now no -- \
    -f $STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk -j 32 \
    > .rna_tophat_make_log.txt \
    2>.rna_tophat_make_err.txt \
    &

