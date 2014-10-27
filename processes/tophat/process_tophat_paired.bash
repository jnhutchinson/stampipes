#FINAL_BAM=${SAMPLE_NAME}.tophat.bam

source $MODULELOAD
module load bowtie/1.0.0
module load java/jdk1.7.0_05
module load bedops/2.4.2
module load tophat/2.0.11
module load picard/1.118
module load cufflinks/2.0.2

export SCRIPT_DIR="$STAMPIPES/scripts/tophat"
export REF_DIR="$STAMPIPES/data/tophat/refseq"

if [ ! -e ${SAMPLE_NAME}_R1_fastqc -o ! -e ${SAMPLE_NAME}_R2_fastqc ]; then
qsub -N ".fq${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/fastqc.mk
__SCRIPT__
fi

# Not ideal, eats up entire node until it's done.
qsub -cwd -V -q all.q -N .th-$SAMPLE_NAME -now no -pe threads 8 <<__MAKE__
    make -f $STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk -j 8
__MAKE__
