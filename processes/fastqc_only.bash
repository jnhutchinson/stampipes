# Dependencies
source $MODULELOAD
module load anaconda/2.1.0-2.7 # Used by data upload

cd $FASTQ_DIR

if [ ! -e ${SAMPLE_NAME}_R1_fastqc -o ! -e ${SAMPLE_NAME}_R2_fastqc ]; then

qsub -N ".fq${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`


  make -f $STAMPIPES/makefiles/fastqc.mk

__SCRIPT__

fi

