#!/usr/bin/env nextflow
/*
 * This is a proof-of-concept workflow for running our DNase alignment pipeline
 */

params.help = false
params.threads = 1
params.chunk_size = 5000  //TODO: raise
params.UMI = false
params.genome = ""

params.outdir = "output"

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
    trim-adapters-illumina \
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
  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${params.threads}" \
    "$genome" \
    "$trimmed_r1" \
    > out1.sai

  bwa aln \
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
  samtools index merged.bam
  """
}

/*
 * Step 4: Mark duplicates with Picard
 */
process mark_duplicates {

  module 'jdk/1.8.0_92'
  module 'picard/2.8.1'
  module 'samtools/1.3'

  memory '40 GB'

  publishDir params.outdir

  input:
  file(merged_bam) from merged_bam

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

/*
 * Step 5: Filter bam file
 */
filter_flag = 512
if (params.UMI)
  filter_flag = 1536

process filter_bam_to_unique {

  module 'samtools/1.3'

  publishDir params.outdir

  input:
  file marked_bam

  output:
  set file('filtered.bam'), file('filtered.bam.bai') into uniquely_mapping_bam

  script:
  """
  samtools view $marked_bam -b -F $filter_flag > filtered.bam
  samtools index filtered.bam
  """

}

// Marked bam file is used by several downstream results
uniquely_mapping_bam.into { bam_for_insert; bam_for_counts; bam_for_spot; bam_for_density }

/*
 * Metrics: Insert size
 */
process insert_size {
  module 'jdk/1.8.0_92'
  module 'picard/2.8.1'
  module 'samtools/1.3'
  module 'R/3.2.5'

  publishDir params.outdir

  input:
  set file(bam), file(bai) from bam_for_insert

  output:
  file 'CollectInsertSizeMetrics.picard'
  file 'CollectInsertSizeMetrics.picard.pdf'

  script:
  """
  samtools idxstats "$bam" \
  | cut -f 1 \
  | grep -v chrM \
  | grep -v chrC \
  | xargs samtools view -b "$bam" \
  > nuclear.bam

  picard CollectInsertSizeMetrics \
    INPUT=nuclear.bam \
    OUTPUT=CollectInsertSizeMetrics.picard \
    HISTOGRAM_FILE=CollectInsertSizeMetrics.picard.pdf \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true
  """
}

/*
 * Metrics: SPOT score
 */

process spot_score {

  module 'samtools/1.3'
  module 'python/2.7.11'
  module 'python/3.5.1'
  module 'pysam/0.9.0'
  module 'bedops/2.4.19'
  module 'bedtools/2.25.0'
  module "R/3.2.5"

  publishDir params.outdir

  input:
  set file(bam), file(bai) from bam_for_spot

  output:
  file 'subsample.spot.out'

  script:
  """
  # random sample
	samtools view -h -F 12 -f 3 "$bam" \
		| awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
		| samtools view -uS - \
		> paired.bam
	bash \$STAMPIPES/scripts/bam/random_sample.sh paired.bam subsample.bam 5000000

  # hotspot
  bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
    /home/solexa/hotspot-hpc/hotspot-distr \
    \$PWD \
    \$PWD/subsample.bam \
    GRCh38_no_alts \
    36 \
    DNaseI
  """
}

win = 75
bini = 20
process density_files {

  module 'bedops/2.4.19'
  module 'samtools/1.3'
  module 'htslib/1.6.0'
  module 'kentutil/302'

  publishDir params.outdir

  input:
  set file(bam), file(bai) from bam_for_density

  output:
  file 'density.bed.starch'
  file 'density.bw'
  file 'density.bed.bgz'


  script:
  """
	bam2bed -d \
	  < $bam \
	  | cut -f1-6 \
    | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > sample.bed

	unstarch "\$STAMPIPES_DATA/densities/chrom-buckets.GRCh38_no_alts.${win}_${bini}.bed.starch" \
    | bedmap --faster --echo --count --delim "\t" - sample.bed \
    | awk -v binI=$bini -v win=$win \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
    | starch - \
    > density.bed.starch

  unstarch density.bed.starch | awk -v binI=$bini -f \$STAMPIPES/awk/bedToWig.awk > density.wig

  wigToBigWig -clip density.wig "${genome}.fai" density.bw


  unstarch density.bed.starch | bgzip > density.bed.bgz
  tabix -p bed density.bed.bgz
  """
}
