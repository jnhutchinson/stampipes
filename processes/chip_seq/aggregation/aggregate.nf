params.help = false
params.threads = 10

params.UMI = false
params.genome = ""
params.control = ""
params.outdir = "output"
params.domotifs = false
params.dofeatures = false

params.readlength = 36

def helpMessage() {
  log.info"""
    Usage: nextflow run basic.nf \\
             --genome /path/to/genome \\
             --treatment '1.bam,2.bam...' \\
             --control 'control.bam(optional)' \\
	     --UMI true/false        \\
             --mappable 'path/to/mappable/file' \\
             --outdir /path/to/output
  """.stripIndent();
}

dataDir = "$baseDir/../../../data"
genome_name = file(params.genome).baseName

bams = Channel.from(
  params.treatment.tokenize(',')
).map {
  file(it)
}.collect()

process merge {

  input:
  file 'in*.bam' from bams

  output:
  file 'merged.bam' into merged

  publishDir params.outdir

  script:
  """
  samtools merge 'merged.bam' in*.bam
  """
}

// TODO: single end
process dups {
  publishDir params.outdir
  memory '16 GB'

  input:
  file(merged)

  output:
  file 'marked.bam' into marked_bam
  file 'MarkDuplicates.picard'

  script:
  if (params.UMI)
    cmd = "UmiAwareMarkDuplicatesWithMateCigar"
    extra = "UMI_TAG_NAME=XD"
  if (!params.UMI)
    cmd = "MarkDuplicatesWithMateCigar"
    extra = ""
  """
  picard RevertOriginalBaseQualitiesAndAddMateCigar \
    "INPUT=${merged}" OUTPUT=cigar.bam \
    VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0
  picard "${cmd}" \
      INPUT=cigar.bam OUTPUT=marked.bam \
      $extra \
      METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
      MINIMUM_DISTANCE=300
  """
}

marked_bam.into { bam_for_counts; bam_for_adapter_counts; bam_for_filter }

process filter {

  publishDir params.outdir

  input:
  file bam from bam_for_filter

  output:
  file "filtered.bam" into filtered_bam

  script:
  flag = params.UMI ? 1536 : 512
  """
  samtools view -b -F "${flag}" marked.bam > filtered.bam
  """
}
filtered_bam.into {  bam_for_macs2; bam_for_spot_score; bam_for_cutcounts; bam_for_inserts; bam_for_nuclear }

process filter_nuclear {
  
  publishDir params.outdir

  input:
  file bam from bam_for_nuclear
  file nuclear_chroms from file("${params.genome}.nuclear.txt")

  output:
  file 'nuclear.bam' into nuclear_bam

  script:
  """
  samtools index "${bam}"
  cat "${nuclear_chroms}" \
  | xargs samtools view -b "${bam}" \
  > nuclear.bam
  """
}

//Control sample to be included?
process macs2_human {
  
  publishDir params.outdir

  input:
  file bam from bam_for_macs2
 
  output:
  file('peaks/filtered*') 
  file('peaks/filtered_treat_pileup.bdg') into signal_treatment
  file('peaks/filtered_control_lambda.bdg') into signal_control

  when:
  genome_name == "GRCh38_no_alts"
     
  script:
  """
    macs2 callpeak -t "${bam}" -f BAM -g hs --outdir peaks -q 0.05 -B -n filtered
  """
}

process macs2_mouse {
  input:
  file bam from bam_for_macs2

  output:
  file('peaks/filtered*')

  when:
  genome_name == "mm10"

  script:
  """
    macs2 callpeak -t "${bam}" -f BAM -g mm --outdir peaks -q 0.05 -B -n filtered
  """
}

process spot_score {
  publishDir params.outdir

  input:
  file(bam) from bam_for_spot_score
  file(mappable) from file(params.mappable)
  //file(mappable) from file("${dataDir}/annotations/${genome_name}.K${params.readlength}.mappable_only.bed")
  file(chromInfo) from file("${dataDir}/annotations/${genome_name}.chromInfo.bed")

  output:
  file 'r1.spot.out'
  file 'r1.hotspot.info'

  script:
  """
  # random sample
        samtools view -h -F 12 -f 3 "$bam" \
                | awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
                | samtools view -uS - \
                > nuclear.bam
        bash \$STAMPIPES/scripts/bam/random_sample.sh nuclear.bam subsample.bam 5000000
  samtools view -b -f 0x0040 subsample.bam > r1.bam
  # hotspot
  bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
    \$HOTSPOT_DIR \
    \$PWD \
    \$PWD/r1.bam \
    "${genome_name}" \
    "${params.readlength}" \
    ChIP_Seq

  starch --header r1-both-passes/r1.hotspot.twopass.zscore.wig \
    > r1.spots.starch

  bash \$STAMPIPES/scripts/SPOT/info.sh \
    r1.spots.starch hotspot1 r1.spot.out \
    > r1.hotspot.info
  """
}

process bam_counts {
  publishDir params.outdir

  input:
  file(bam) from bam_for_counts

  output:
  file('tagcounts.txt')

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/bamcounts.py \
    "$bam" \
    tagcounts.txt
  """
}

process count_adapters {
  publishDir params.outdir

  input:
  file(bam) from bam_for_adapter_counts

  output:
  file('adapter.counts.txt')

  script:
  """
  bash "\$STAMPIPES/scripts/bam/count_adapters.sh" "${bam}" \
  | sed 's/^/adapter\t/' \
  > adapter.counts.txt
  """
}

process preseq {
  publishDir params.outdir
  input:
  file nuclear_bam

  output:
  file 'preseq.txt'
  file 'preseq_targets.txt'
  file 'dups.hist'

  script:
  """
  python3 \$STAMPIPES/scripts/bam/mark_dups.py -i "${nuclear_bam}" -o /dev/null --hist dups.hist
  preseq lc_extrap -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt \
  || preseq lc_extrap -defects -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt

  # write out preseq targets
  bash "\$STAMPIPES/scripts/utility/preseq_targets.sh" preseq.txt preseq_targets.txt
  """
}



process signal_compare {
   publishDir params.outdir
  
   input : 
   file treatment from signal_treatment
   file control from signal_control

   output : 
   file('peaks/signal_qpois.bw')  

   script : 
   """
   macs2 bdgcmp -t "${treatment}" -c "${control}" --outdir peaks --o-prefix signal -m qpois 
   bedtools slop -i peaks/signal_qpois.bdg -g "${dataDir}/annotations/GRCh38_no_alts.chrom.sizes" -b 0 | bedClip stdin "${dataDir}/annotations/GRCh38_no_alts.chrom.sizes" peaks/signal_qpois.bdg.clip
   LC_COLLATE=C sort -k1,1 -k2,2n peaks/signal_qpois.bdg.clip > peaks/signal_qpois.bdg.clip.sort
   bedGraphToBigWig peaks/signal_qpois.bdg.clip.sort "${dataDir}/annotations/GRCh38_no_alts.chrom.sizes" peaks/signal_qpois.bw
   rm peaks/signal_qpois.bdg.clip.sort peaks/signal_qpois.bdg.clip
   """
}

process insert_sizes {

  publishDir params.outdir

  input:
  file nuclear_bam from bam_for_inserts
  file nuclear_chroms from file("${params.genome}.nuclear.txt")

  output:
  file 'CollectInsertSizeMetrics.picard*'

  script:
  """
  picard CollectInsertSizeMetrics \
    "INPUT=${nuclear_bam}" \
    OUTPUT=CollectInsertSizeMetrics.picard \
    HISTOGRAM_FILE=CollectInsertSizeMetrics.picard.pdf \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true

  cat CollectInsertSizeMetrics.picard \
  | awk '/## HISTOGRAM/{x=1;next}x' \
  | sed 1d \
  > hist.txt

  python3 "\$STAMPIPES/scripts/utility/picard_inserts_process.py" hist.txt > CollectInsertSizeMetrics.picard.info
  """
}
