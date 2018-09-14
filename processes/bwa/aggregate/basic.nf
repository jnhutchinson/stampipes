params.help = false
params.threads = 1

params.UMI = false
params.genome = ""
params.outdir = "output"
params.domotifs = false
params.dofeatures = false

params.readlength = 36

def helpMessage() {
  log.info"""
    Usage: nextflow run basic.nf \\
             --genome /path/to/genome \\
             --bams '1.bam,2.bam...' \\
             --UMI true/false        \\
             --outdir /path/to/output

  """.stripIndent();
}

dataDir = "$baseDir/../../../data"
genome_name = file(params.genome).baseName

bams = Channel.from(
  params.bams.tokenize(',')
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
filtered_bam.into { bam_for_hotspot2; bam_for_spot_score; bam_for_cutcounts; bam_for_density; bam_for_inserts; bam_for_nuclear }

process filter_nuclear {
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

process hotspot2 {

  publishDir "${params.outdir}"
  container "fwip/hotspot2:latest"

  input:
  file(marked_bam) from bam_for_hotspot2
  file(mappable) from file(params.mappable)
  file(chrom_sizes) from file(params.chrom_sizes)
  file(centers) from file(params.centers)


  output:
  file('peaks/filtered*')
  file('peaks/filtered.hotspots.fdr0.05.starch') into hotspot_calls

  script:
  """
  export TMPDIR=\$(mktemp -d)
  hotspot2.sh -F 0.5 -p varWidth_20_default \
    -M "${mappable}" \
    -c "${chrom_sizes}" \
    -C "${centers}" \
    "${marked_bam}" \
    'peaks'

  cd peaks

  # Rename peaks files to include FDR
  mv filtered.peaks.narrowpeaks.starch filtered.peaks.narrowpeaks.fdr0.05.starch
  mv filtered.peaks.starch filtered.peaks.fdr0.05.starch

  # TODO: Move this to separate process
  hsmerge.sh -f 0.01 filtered.allcalls.starch filtered.hotspots.fdr0.01.starch
  hsmerge.sh -f 0.001 filtered.allcalls.starch filtered.hotspots.fdr0.001.starch

  density-peaks.bash \$TMPDIR varWidth_20_default filtered.cutcounts.starch filtered.hotspots.fdr0.01.starch ../"${chrom_sizes}" filtered.density.starch filtered.peaks.fdr0.01.starch \$(cat filtered.cleavage.total)
  density-peaks.bash \$TMPDIR varWidth_20_default filtered.cutcounts.starch filtered.hotspots.fdr0.001.starch ../"${chrom_sizes}" filtered.density.starch filtered.peaks.fdr0.001.starch \$(cat filtered.cleavage.total)

  rm -rf "\$TMPDIR"
  """

}

process spot_score {
  publishDir params.outdir

  input:
  file(bam) from bam_for_spot_score
  file(mappable) from file("${dataDir}/annotations/${genome_name}.K${params.readlength}.mappable_only.bed")
  file(chromInfo) from file("${dataDir}/annotations/${genome_name}.chromInfo.bed")

  output:
  file 'r1.spot.out'

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
    DNaseI
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
  bash "\$STAMPIPES/scripts/bam/count_adapters.sh" "${bam}" > adapter.counts.txt
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
  preseq lc_extrap -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt

  # write out preseq targets
  bash "\$STAMPIPES/scripts/utility/preseq_targets.sh" preseq.txt preseq_targets.txt
  """
}

process cutcounts {

  publishDir params.outdir

  input:
  file(fai) from file("${params.genome}.fai")
  file(filtered_bam) from bam_for_cutcounts

  output:
  file('fragments.starch')
  file('cutcounts.starch')
  file('cutcounts.bw')
  file('cutcounts.bed.bgz')
  file('cutcounts.bed.bgz.tbi')

  script:
  """
  bam2bed --do-not-sort \
  < "$filtered_bam" \
  | awk -v cutfile=cuts.bed -v fragmentfile=fragments.bed -f \$STAMPIPES/scripts/bwa/aggregate/basic/cutfragments.awk

  sort-bed fragments.bed | starch - > fragments.starch
  sort-bed cuts.bed | starch - > cuts.starch

  unstarch cuts.starch \
  | cut -f1-3 \
  | bedops -m - \
  | awk '{ for(i = \$2; i < \$3; i += 1) { print \$1"\t"i"\t"i + 1 }}' \
  > allbase.tmp

  unstarch cuts.starch \
  | bedmap --echo --count --delim "\t" allbase.tmp - \
  | awk '{print \$1"\t"\$2"\t"\$3"\tid-"NR"\t"\$4}' \
  > cutcounts.tmp

  starch cutcounts.tmp > cutcounts.starch

  awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,\$5}' cutcounts.tmp \
  | grep -v chrM \
  > wig.tmp

  wigToBigWig -clip wig.tmp "${fai}" cutcounts.bw

  # tabix
  unstarch cutcounts.starch | bgzip > cutcounts.bed.bgz
  tabix -p bed cutcounts.bed.bgz
  """
}

process density {

  publishDir params.outdir

  input:
  file filtered_bam from bam_for_density
  file chrom_bucket from file(params.chrom_bucket)

  output:
  file 'density.starch'

  script:
  """
  mkfifo density.bed

  bam2bed -d \
  < "${filtered_bam}" \
  | cut -f1-6 \
  | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
  | sort-bed - \
  > density.bed \
  &

  unstarch "${chrom_bucket}" \
  | bedmap --faster --echo --count --delim "\t" - density.bed \
  | awk -v binI=20 -v win=75 \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
  | starch - \
  > density.starch
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

process motif_matrix {

  publishDir params.outdir

  input:
  file hotspot_calls
  file fimo_transfac from file("${dataDir}/motifs/${genome_name}.fimo.starch")
  file fimo_names from file("${dataDir}/motifs/${genome_name}.fimo.transfac.names.txt")

  output:
  file 'hs_motifs*.txt'

  when:
  params.domotifs

  script:
  """
  # create sparse motifs
  bedmap --echo --echo-map-id --fraction-map 1 --delim '\t' "${hotspot_calls}" "${fimo_transfac}" > temp.bedmap.txt
  python \$STAMPIPES/scripts/bwa/aggregate/basic/sparse_motifs.py "${fimo_names}" temp.bedmap.txt
  """
}

process closest_features {

  publishDir params.outdir

  input:
  file hotspot_calls
  file transcript_starts from file("${dataDir}/features/${genome_name}.CombinedTxStarts.bed")
  val thresholds from "0 1000 2500 5000 10000"

  when:
  params.dofeatures

  output:
  file 'prox_dist.info'

  script:
  """
  closest-features \
    --dist \
    --delim '\t' \
    --closest \
    "${hotspot_calls}" \
    "${transcript_starts}" \
  > closest.txt
  cat closest.txt \
  | grep -v "NA\$" \
  | awk -F"\t" '{print \$NF}' \
  | sed -e 's/-//g' \
  > closest.clean.txt

  for t in $thresholds ; do
    awk \
      -v t=\$t \
      '\$1 > t {sum+=1} END {print "percent-proximal-" t "bp " sum/NR}' \
      closest.clean.txt \
    >> prox_dist.info
  done
  """
}
