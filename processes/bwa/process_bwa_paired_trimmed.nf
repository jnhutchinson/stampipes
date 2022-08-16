nextflow.enable.dsl=2

/*
 * This is our main DNase alignment pipeline
 *
 * It performs several functions:
 * 1) Pre-processing reads (e.g: removing adapters)
 * 2) Alignment to a reference genome with BWA
 * 3) Duplicate marking and removal of low-quality reads
 * 4) QC stat calculation
 */

params.help = false
params.threads = 1
params.chunk_size = 16000000

params.id = ""
params.UMI = false
params.trim_to = 0
params.genome = ""
params.r1 = ""
params.r2 = ""
params.outdir = "output"
params.readlength = 36
params.adapter_file = ""
//not sure if this would work with single end data with paired=true as default
params.paired = false

params.cramthreads = 10

def helpMessage() {
  log.info"""
    Usage: nextflow run process_bwa_paired_trimmed.nf \\
             --r1 r1.fastq.gz \\
             --r2 r2.fastq.gz \\
             --adapter_file adapters.txt \\
             --genome /path/to/genome \\
    Options:
    --threads [count]     The number of threads that will be used for child applications  (1)
    --chunk_size [count]  How many reads to process per chunk                             (16000000)
    --UMI                 The reads contain UMI markers ('single-strand', 'thruplex')         (false)
    --trim_to [length]    Trim fastq reads to [length] bp (0 for no trimming)             (0)
    --outdir [dir]        Where to write results to                                       (output)
    """.stripIndent();
}

// Used for publishDir to avoid publishing bam files
def exclude_bam =  { name -> name ==~ /.*ba[mi]$/ ? null : name }


if (params.help || !params.r1 ||  !params.genome){
  helpMessage();
  exit 0;
}

// Some renaming for easier usage later
dataDir = "${baseDir}/../../data"

include { split_paired } from "../../modules/fastq.nf" addParams(fastq_split_count: params.chunk_size)
//include { get_bwa_index_files; bwa_sampe } from "../../modules/bwa.nf"
include { BWA_ALIGN } from "../../modules/bwa.nf"
include { fastp_adapter_trim; parse_legacy_adapter_file } from "../../modules/adapter_trimming.nf"
include { encode_cram } from "../../modules/cram.nf"
include { publish; publish_with_meta } from "../../modules/utility.nf"

// When this file is called directly, run our default BWA_ALIGNMENT pipeline for a single sample.
workflow {
  def metadata = [
    id: params.id,
    r1: file(params.r1, checkIfExists: true),
    r2: file(params.r2),
    genome: file(params.genome, checkIfExists: true),
    adapter_file: file(params.adapter_file, checkIfExists: true),
    read_length: params.readlength,
    umi: params.UMI,
    paired: params.paired,
    outdir: params.outdir,
  ]
  BWA_ALIGNMENT(channel.of(metadata))
}

// Fill out default values if not provided to us
def meta_defaults(m) {
  def meta = m.clone()
  meta.putIfAbsent("trim_to_length", 0)
  meta.putIfAbsent("genome_name", meta.genome.simpleName)
  meta.putIfAbsent("density_window_width", 75)
  meta.putIfAbsent("density_step_size", 20)

  meta.putIfAbsent("reference_nuclear_chroms", file("${meta.genome}.nuclear.txt", checkIfExists: true))
  meta.putIfAbsent("reference_fai", file("${meta.genome}.fai", checkIfExists: true))
  meta.putIfAbsent("reference_chrom_info", file(
    "${dataDir}/annotations/${meta.genome_name}.chromInfo.bed",
    checkIfExists: true
  ))
  meta.putIfAbsent("reference_chrom_buckets", file(
    "${dataDir}/densities/chrom-buckets.${meta.genome_name}.${meta.density_window_width}_${meta.density_step_size}.bed.starch",
    checkIfExists: true
  ))
  meta.putIfAbsent("reference_mappable", file(
    "${dataDir}/annotations/${meta.genome_name}.K${meta.read_length}.mappable_only.bed",
    checkIfExists: true
  ))
  return meta
}


workflow BWA_ALIGNMENT {
  take: orig_metadata

  main:
    def metadata = orig_metadata.map { meta_defaults(it) }

    metadata.map { meta -> [meta, file(meta.r1)] }
    | fastq_counts

    metadata.map { meta -> [meta, file(meta.r1), file(meta.r2)] }
    //| view {it, r1, r2 -> "input to split_paired: ${it}, ${r1}, ${r2}"}
    | split_paired
    //| view { "to map $it" }
    | map { m, reads ->
         (p7, p5) = parse_legacy_adapter_file(m.adapter_file);
         return [m, reads[0], reads[1], p7, p5]
       }
    // | view { "from map $it" }
    | fastp_adapter_trim

  //if single end data, submit null R2_fastq to BWA for alignment
    if (params.paired){   
    fastp_adapter_trim.out.fastq
      .map { meta, reads -> [ meta, meta.trim_to_length, reads[0], reads[1] ] }
    | trim_to_length
    | map { meta, r1, r2 -> [meta, meta.umi, r1, r2] }
    | add_umi_info
    | map { meta, r1, r2 -> [ meta, r1, r2, meta.genome ]  }
    | BWA_ALIGN
    | map { meta, bam -> [ meta, bam, meta.reference_nuclear_chroms ] }
    | filter_bam
    | sort_bam
    | groupTuple()  // TODO: This blocks until all sort_bams are done... seems bad.
    | merge_bam
    | map { meta, bam -> [ meta, !!meta.umi, bam ] }
    | mark_duplicates
    } else {
     fastp_adapter_trim.out.fastq
      .map { meta, reads -> [ meta, meta.trim_to_length, reads[0], reads[1] ] }
     | trim_to_length
     | map { meta, r1, r2 -> [meta, meta.umi, r1, r2] }
     | add_umi_info
     | map { meta, r1, r2 -> [ meta, r1, "", meta.genome ]  }
     | BWA_ALIGN
     | map { meta, bam -> [ meta, bam, meta.reference_nuclear_chroms ] }
     | filter_bam
     | sort_bam
     | groupTuple()  // TODO: This blocks until all sort_bams are done... seems bad.
     | merge_bam
     | map { meta, bam -> [ meta, !!meta.umi, bam ] }
     | mark_duplicates
    }

    mark_duplicates.out.marked_bam
      .map { meta, bam -> [ meta, !!meta.umi, bam ] }
    | filter_bam_to_unique

    mark_duplicates.out.marked_bam
    | bam_counts

  //skip insert size calculation if single end data
   if (params.paired){   
    filter_bam_to_unique.out
    | insert_size
   }

    filter_bam_to_unique.out
      .map { meta, bam, bai -> [ meta, bam, bai,
        meta.reference_mappable,
        meta.reference_chrom_info,
      ]}
    | view { "b4 spot $it" }
    | spot_score

    filter_bam_to_unique
      .out.map { meta, bam, bai -> [ meta, bam, bai,
        meta.reference_fai,
        meta.reference_chrom_buckets,
      ]}
    | density_files

    //filter_bam_to_unique
    //  .out.map { meta, bam, bai -> [meta, bam, meta.genome] }
    //| encode_cram

    mark_duplicates.out.marked_bam
    //filter_bam_to_unique
      .map { meta, bam -> [meta, bam, meta.genome] }
    | encode_cram

    Channel.empty().mix(
      // fastp_adapter_trim.out.counts,
      bam_counts.out,
      fastq_counts.out,
    ).groupTuple()
    | gather_counts

    //skip publishing insert size if single end 
   if (params.paired){   
    // Publish all files
    Channel.empty().mix(
      mark_duplicates.out.metrics,
      insert_size.out.plaintext,
      gather_counts.out,
      density_files.out.starch,
      density_files.out.bigwig,
      density_files.out.bgzip,
    ) | publish_with_meta
   } else {
     // Publish all files
    Channel.empty().mix(
      mark_duplicates.out.metrics,
      gather_counts.out,
      density_files.out.starch,
      density_files.out.bigwig,
      density_files.out.bgzip,
    ) | publish_with_meta 
   }

    Channel.empty().mix(
      encode_cram.out.cram_with_index.flatMap { meta, cram, crai -> [cram, crai] },
      spot_score.out.flatMap { meta, score, dups -> [score, dups] }
    ) | publish

}

process sort_bam {
  cpus params.threads

  input:
    tuple val(meta), file(filtered_bam)

  output:
    tuple val(meta), file('sorted.bam')

  script:
  """
  samtools sort \
    -l 1 -m 1G -@ "${params.threads}" "${filtered_bam}" \
    > sorted.bam
  """
}

process add_umi_info {

  input:
  tuple val(meta), val(UMI), path(r1), path(r2)

  output:
  tuple val(meta), file('r1.fastq.umi.gz'), file('r2.fastq.umi.gz')

  script:
  if (UMI == 'thruplex')
    """
    python3 "\$STAMPIPES/scripts/umi/extract_umt.py" \
      <(zcat "${r1}") \
      <(zcat "${r2}") \
      >(gzip -c -1 > r1.fastq.umi.gz) \
      >(gzip -c -1 > r2.fastq.umi.gz)
    """
  else if (UMI == 'single-strand')
    """
    python3 "\$STAMPIPES/scripts/umi/fastq_umi_add.py" "${r1}" r1.fastq.umi.gz
    python3 "\$STAMPIPES/scripts/umi/fastq_umi_add.py" "${r2}" r2.fastq.umi.gz
    """

  else if (UMI == false || UMI == "" || UMI == null)
    """
    ln -s "${r1}" r1.fastq.umi.gz
    ln -s "${r2}" r2.fastq.umi.gz
    """
  else
    error "--UMI must be `thruplex`, `single-strand` (for single-strand preparation), or false, got: '${UMI}'"
}

process trim_to_length {

  input:
  tuple val(meta), val(length), path(r1), path(r2)

  output:
  tuple val(meta), file('r1.trim.fastq.gz'), file('r2.trim.fastq.gz')

  script:
  if (length > 0)
    // TODO: Add padding to length with N's
    """
    zcat "${r1}" | awk 'NR%2==0 {print substr(\$0, 1, ${length})} NR%2!=0' | gzip -c -1 > r1.trim.fastq.gz
    zcat "${r2}" | awk 'NR%2==0 {print substr(\$0, 1, ${length})} NR%2!=0' | gzip -c -1 > r2.trim.fastq.gz
    """
  else
    """
    ln -s "${r1}" r1.trim.fastq.gz
    ln -s "${r2}" r2.trim.fastq.gz
    """

}

process merge_bam {
  input:
    tuple val(meta), path('sorted_bam_*')

  output:
    tuple val(meta), file('merged.bam')

  script:
    """
    samtools merge merged.bam sorted_bam_*
    """
}

process mark_duplicates {

  label "high_mem"

  input:
    tuple val(meta), val(umi_aware), path(merged_bam)

  output:
    tuple val(meta), file('marked.bam'), emit: marked_bam
    tuple val(meta), file('MarkDuplicates.picard'), emit: metrics


  script:
    // TODO: Why do we use MINIMUM_DISTANCE only if we don't have UMIs?
    if (umi_aware)
      """
      picard RevertOriginalBaseQualitiesAndAddMateCigar \
        INPUT="${merged_bam}" OUTPUT=cigar.bam \
        VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

      picard UmiAwareMarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
        METRICS_FILE=MarkDuplicates.picard UMI_TAG_NAME=XD ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
      """
    else
      """
      picard RevertOriginalBaseQualitiesAndAddMateCigar \
        INPUT="${merged_bam}" OUTPUT=cigar.bam \
        VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

      picard MarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
        METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
        MINIMUM_DISTANCE=300
      """
}

process filter_bam_to_unique {
  input:
    tuple val(meta), val(is_umi), path(marked_bam)

  output:
    tuple val(meta), file('filtered.bam'), file('filtered.bam.bai')

  script:
    filter_flag = is_umi ? 1536 : 512
    """
    samtools view "${marked_bam}" -b -F "${filter_flag}" -o filtered.bam
    samtools index filtered.bam
    """
}

process insert_size {
  // TODO: Fix chomosome removal!!
  input:
    tuple val(meta), path(bam), path(bai)

  output:
    tuple val(meta), file('CollectInsertSizeMetrics.picard'), emit: plaintext
    tuple val(meta), file('CollectInsertSizeMetrics.picard.pdf'), emit: pdf

  script:
    """
    samtools idxstats "${bam}" \
    | cut -f 1 \
    | grep -v chrM \
    | grep -v chrC \
    | xargs samtools view -b "${bam}" -o nuclear.bam

    picard CollectInsertSizeMetrics \
      INPUT=nuclear.bam \
      OUTPUT=CollectInsertSizeMetrics.picard \
      HISTOGRAM_FILE=CollectInsertSizeMetrics.picard.pdf \
      VALIDATION_STRINGENCY=LENIENT \
      ASSUME_SORTED=true
    """
}

process bam_counts {
  input:
    tuple val(meta), file(sorted_bam)

  output:
    tuple val(meta), file('bam.counts.txt')

  script:
    """
    python3 \$STAMPIPES/scripts/bwa/bamcounts.py \
      "$sorted_bam" \
      bam.counts.txt
    """
}

process spot_score {
 scratch false
  input:
    tuple val(meta), path(bam), path(bai), path(mappable_only), path(chrom_info)

  output:
    tuple val(meta), file('subsample.r1.spot.out'), file('spotdups.txt')

  script:
    def genome_name = chrom_info.simpleName
    def read_length = (mappable_only.name =~ /K([0-9]+)/)[0][1]
  // conditional to deal with single end data
  if (params.paired)  
    """
    # random sample
    samtools view -h -F 12 -f 3 "${bam}" \
      | awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
      | samtools view -1 - \
      -o paired.bam
    bash \$STAMPIPES/scripts/bam/random_sample.sh paired.bam subsample.bam 5000000
          samtools view -1 -f 0x0040 subsample.bam -o subsample.r1.bam

    # hotspot
    bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
      \$HOTSPOT_DIR \
      \$PWD \
      \$PWD/subsample.r1.bam \
      "${genome_name}" \
      "${read_length}" \
      DNaseI

    # Remove existing duplication marks
    picard RevertSam \
      INPUT=subsample.bam \
      OUTPUT=clear.bam \
      VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATE_INFORMATION=true SORT_ORDER=coordinate \
      RESTORE_ORIGINAL_QUALITIES=false REMOVE_ALIGNMENT_INFORMATION=false

    picard MarkDuplicatesWithMateCigar \
      INPUT=clear.bam \
      METRICS_FILE=spotdups.txt \
      OUTPUT=/dev/null \
      ASSUME_SORTED=true \
      MINIMUM_DISTANCE=300 \
      VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """
    else
    """
    # random sample
    samtools view -h -F 12 -s 0.6 "${bam}" \
      | awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
      | samtools view -1 - \
      -o paired.bam
    bash \$STAMPIPES/scripts/bam/random_sample.sh paired.bam subsample.bam 5000000
          samtools view -1 subsample.bam -o subsample.r1.bam

    # hotspot
    bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
      \$HOTSPOT_DIR \
      \$PWD \
      \$PWD/subsample.r1.bam \
      "${genome_name}" \
      "${read_length}" \
      DNaseI
    
    # Remove existing duplication marks
    picard RevertSam \
      INPUT=subsample.bam \
      OUTPUT=clear.bam \
      VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATE_INFORMATION=true SORT_ORDER=coordinate \
      RESTORE_ORIGINAL_QUALITIES=false REMOVE_ALIGNMENT_INFORMATION=false

    picard MarkDuplicates \
      INPUT=clear.bam \
      METRICS_FILE=spotdups.txt \
      OUTPUT=/dev/null \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' 
    """
}

process density_files {

  label "high_mem"

  input:
    tuple val(meta), path(bam), path(bai), path(fai), path(density_buckets)

  output:
    tuple val(meta), file('density.bed.starch'), emit: starch
    tuple val(meta), file('density.bw'), emit: bigwig
    tuple val(meta), file('density.bed.bgz'), emit: bgzip


  script:
    win = meta.density_window_width
    bini = meta.density_step_size
    """
    bam2bed -d \
      < "${bam}" \
      | cut -f1-6 \
      | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
      | sort-bed - \
      > sample.bed

    unstarch "${density_buckets}" \
      | bedmap --faster --echo --count --delim "\t" - sample.bed \
      | awk -v binI=$bini -v win="${win}" \
          'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
      | starch - \
      > density.bed.starch

    unstarch density.bed.starch | awk -v binI="${bini}" -f "\$STAMPIPES/awk/bedToWig.awk" > density.wig

    wigToBigWig -clip density.wig "${fai}" density.bw

    unstarch density.bed.starch | bgzip > density.bed.bgz
    tabix -p bed density.bed.bgz
    """
}

process gather_counts {
  input:
    tuple val(meta), file("counts*.txt")

  output:
    tuple val(meta), file('tagcounts.txt')

  script:
    """
    cat counts*.txt \
    | awk '
      { x[\$1] += \$2 }
      END {for (i in x) print i "\t" x[i]}
    ' \
    | sort -k 1,1 \
    > tagcounts.txt
    """
}

process filter_bam {

  input:
  tuple val(meta), path(unfiltered_bam), path(nuclear_chroms)

  output:
  tuple val(meta), file('filtered.bam')

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/filter_reads.py \
  "${unfiltered_bam}" \
  filtered.bam \
  "${nuclear_chroms}"
  """
}

process fastq_counts {

  input:
    tuple val(meta), path(r1)

  output:
    tuple val(meta), file('fastq.counts')

  script:
    """
    zcat "${r1}" \
    | awk -v paired=1 -f \$STAMPIPES/awk/illumina_fastq_count.awk \
    > fastq.counts
    """
}
