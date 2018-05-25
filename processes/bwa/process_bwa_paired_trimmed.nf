#!/usr/bin/env nextflow
/*
 * This is a proof-of-concept workflow for running our DNase alignment pipeline
 */

params.help = false
params.threads = 1
params.chunk_size = 5000  //TODO: raise
params.UMI = false
params.genome = ""

nuclear_chroms = "$params.genome" + ".nuclear.txt"

def helpMessage() {
  log.info"""
    Usage: nextflow run process_bwa_paired_trimmed.nf --r1 r1.fastq.gz --r2 r2.fastq.gz --adapter_file adapters.txt --genome /path/to/genome
    Options:
    --threads [count]     The number of threads that will be used for child applications  (1)
    --chunk_size [count]  How many reads to process per chunk                             (5000)
    --UMI                 The reads contain UMI markers                                   (false)
    """.stripIndent();
}

if (params.help){
  helpMessage();
  exit 0;
}

// Split input files into chunks to process in parallel
split_r1 = Channel
  .from( file(params.r1) )
  .splitFastq( by: params.chunk_size)
split_r2 = Channel
  .from( file(params.r2) )
  .splitFastq( by: params.chunk_size)


// Some renaming for easier usage later
genome = params.genome
threads = params.threads
adapters = file(params.adapter_file)


/*
 * Step 1: For each fastqc chunk, trim adapters
 */
process trim_adapters {

    input:
    file split_r1
    file split_r2

    output:
    file 'trim.R1.fastq.gz' into trimmed_r1
    file 'trim.R2.fastq.gz' into trimmed_r2

    script:
    """
    time trim-adapters-illumina \
      -f "$adapters" \
      -1 P7 -2 P5 \
      --threads=${params.threads} \
      "$split_r1" \
      "$split_r2"  \
      "trim.R1.fastq.gz" \
      "trim.R2.fastq.gz"
    """
}

/*
 * Step 2a: Create alignment files
 */
process align {

  module 'bwa/0.7.12'
  module 'samtools/1.3'

  input:
  file trimmed_r1
  file trimmed_r2

  output:
  file 'out.bam' into unfiltered_bam

  script:
  """
  time bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${params.threads}" \
    "$genome" \
    "$trimmed_r1" \
    > out1.sai

  time bwa aln \
    -Y -l 32 -n 0.04 \
    -t "$params.threads" \
    "$genome" \
    "$trimmed_r2" \
    > out2.sai

  bwa sampe \
    -n 10 -a 750 \
    "$genome" \
    out1.sai out2.sai \
    "$trimmed_r1" "$trimmed_r2" \
  | samtools view -b -t "$genome".fai - \
  > out.bam
  """

}

/*
 * Step 2b: filter bam files to have only good reads
 */
process filter_bam {
  module 'samtools/1.3';
  module 'python/3.5.1';
  module 'pysam/0.9.0';

  input:
  file unfiltered_bam

  output:
  file 'filtered.bam' into filtered_bam

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/filter_reads.py \
  "$unfiltered_bam" \
  filtered.bam \
  "$nuclear_chroms"
  """
}

/*
 * Step 2c: sort bam files
 */
process sort_bam {
  module 'samtools/1.3';

  input:
  file filtered_bam

  output:
  file 'sorted.bam' into sorted_bam

  script:
  """
  samtools sort \
    -l 0 -m 1G -@ "${params.threads}" "$filtered_bam" \
    > sorted.bam
  """
}

/*
 * Step 3: Merge alignments into one big ol' file
 */
process merge_bam {
  module 'samtools/1.3'

  input:
  file 'sorted_bam_*' from sorted_bam.collect()

  output:
  file 'merged.bam' into merged_bam

  script:
  """
  samtools merge merged.bam sorted_bam*
  """
}

/*
 * Step 4: Mark duplicates with Picard
 */
process mark_duplicates {

  module 'jdk/1.8.0_92'
  module 'picard/2.8.1'

  memory '40 GB'

  publishDir 'dup_marked'


  input:
  file merged_bam

  output:
  file 'marked.bam' into marked_bam
  file 'metrics.duplicate.picard' into duplicate_report


  script:
  if (params.UMI)
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      INPUT=$merged_bam OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    picard UmiAwareMarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
      METRICS_FILE=metrics.duplicate.picard UMI_TAG_NAME=XD ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """
  else
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      INPUT=$merged_bam OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    picard MarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
      METRICS_FILE=metrics.duplicate.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """
}
