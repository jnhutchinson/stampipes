# Dependencies
source $MODULELOAD
module load java/jdk1.7.0_05
module load picard/1.118
module load FastQC/0.11.3

cd $FASTQ_DIR

if [ ! -e ${SAMPLE_NAME}_R1_fastqc.zip -o ! -e ${SAMPLE_NAME}_R2_fastqc.zip ]; then

PROCESS_NAME=".fq${SAMPLE_NAME}_${FLOWCELL}_LANE#${FLOWCELL_LANE_ID}"

qsub -N ".fq${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  cd $FASTQ_DIR
  make -f $STAMPIPES/makefiles/fastqc.mk

  if [ "$UMI" = "True" ]; then
      echo "Tallying up top UMI tags seen in R1"
      zcat ${R1_FASTQ} | grep "^@" | cut -f 2 -d "+" | sort | uniq -c | sort -n -r | head -n 100 > ${SAMPLE_NAME}.topumis.txt
  fi

  bash $STAMPIPES/scripts/fastq/attachfiles.bash

  echo "FINISH: "
  date
__SCRIPT__
fi

