nextflow.enable.dsl = 2

params.help = false
params.threads = 1

params.UMI = false
params.genome = ""
params.outdir = "output"
params.domotifs = false
params.dofeatures = false

params.readlength = 36

params.peakcaller = "hotspot2"

params.bias = ""
params.chunksize = 5000

params.hotspot_id = "default"
params.hotspot_index = "."
params.bams = ""

params.cramthreads = 10

params.mappable = ""
params.centers = ""
params.chrom_sizes = ""

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

include { encode_cram } from "../../../modules/cram.nf"
include { publish; publish_with_meta; publish_many } from "../../../modules/utility.nf"

workflow {
  def bams = params.bams.tokenize(",").collect { file(it, checkExists: true) }.flatten()
  def meta = [
    genome: file(params.genome, checkExists: true),
    bams: bams,
    UMI:  params.UMI,
    read_length: params.readlength,

    id: params.id,
    hotspot_id: params.hotspot_id ?: params.id,
    hotspot_index: params.hotspot_index,

    reference_mappable: file(params.mappable, checkExists: true),
    reference_chrom_sizes: file(params.chrom_sizes, checkExists: true),
    reference_hotspot_centers: file(params.centers, checkExists: true),

    outdir: params.outdir,
  ]

  DNASE_AGGREGATION(channel.of(meta))
}

def meta_defaults(meta) {
  m = meta.clone()

  // Put new stuff here
  m.putIfAbsent("UMI", false)
  m.putIfAbsent("reference_fai", file("${m.genome}.fai"))
  m.putIfAbsent("genome_name", m.genome.simpleName)
  m.putIfAbsent("density_window_width", 75)
  m.putIfAbsent("density_step_size", 20)
  m.putIfAbsent("reference_nuclear_chroms", file("${m.genome}.nuclear.txt", checkIfExists: true))
  m.putIfAbsent("reference_chrom_buckets", file(
    "${dataDir}/densities/chrom-buckets.${m.genome_name}.${m.density_window_width}_${m.density_step_size}.bed.starch",
    checkIfExists: true
  ))
  m.putIfAbsent("reference_mappable", file(
    "${dataDir}/annotations/${m.genome_name}.K${m.read_length}.mappable_only.bed",
    checkIfExists: true
  ))
  m.putIfAbsent("reference_chrom_info", file(
    "${dataDir}/annotations/${m.genome_name}.chromInfo.bed", checkIfExists: true
  ))

  m.putIfAbsent("fimo_names", file(
    "${dataDir}/motifs/${m.genome_name}.fimo.transfac.names.txt",
    checkIfExists: true))
  m.putIfAbsent("fimo_transfac", file(
    "${dataDir}/motifs/${m.genome_name}.fimo.starch",
    checkIfExists: true))
  m.putIfAbsent("transcript_starts", file(
    "${dataDir}/features/${m.genome_name}.CombinedTxStarts.bed",
    checkIfExists: true))
  return m
}

workflow DNASE_AGGREGATION {

  take: orig_metadata

  main:
    def metadata = orig_metadata.map { meta_defaults(it) }

    metadata.map { meta -> [meta, meta.bams] }
    | merge_bam

    merge_bam.out.map { meta, bam -> [meta, meta.UMI, bam] }
    | dups

    dups.out.bam | bam_counts
    dups.out.bam | count_adapters
    dups.out.bam | map { meta, bam -> [meta, meta.reference_chrom_buckets, meta.reference_fai, bam] }
    | multimapping_density
    dups.out.bam | filter_bam

    filter_bam.out.map { meta, bam -> [meta, meta.reference_chrom_buckets, meta.reference_fai, bam] }
    | density
    filter_bam.out.map { meta, bam -> [meta, meta.reference_fai, bam] }
    | cutcounts
    filter_bam.out.map { meta, bam -> [meta, meta.reference_nuclear_chroms, bam] }
    | insert_sizes
    filter_bam.out.map { meta, bam -> [meta, meta.reference_nuclear_chroms, bam] }
    | filter_nuclear
    filter_bam.out.map { meta, bam -> [meta, bam, meta.genome] }
    | encode_cram

     density.out.starch
       .join( filter_bam.out )
       .map { meta, density, bam -> [meta, meta.reference_fai, density, bam] }
     | normalize_density
 
     filter_nuclear.out | macs2
     filter_nuclear.out | preseq
     filter_nuclear.out.map { meta, bam -> [
       meta, meta.genome_name, meta.reference_mappable, meta.reference_chrom_info, bam] }
     | spot_score
     filter_nuclear.out.map { meta, bam -> [
       meta, meta.hotspot_id, meta.reference_mappable, meta.reference_chrom_sizes, meta.reference_hotspot_centers, bam]} 
     | hotspot2

    hotspot2.out.hotspots \
      | map { meta, hotspots -> [meta, meta.fimo_transfac, meta.fimo_names, hotspots] }
      | motif_matrix
    hotspot2.out.hotspots \
      | map { meta, hotspots -> [meta, meta.transcript_starts, hotspots] }
      | closest_features
    hotspot2.out.peaks
      | join(dups.out.bam)
      | map { meta, peaks, bam -> [meta, meta.hotspot_index, bam, peaks] }
      | differential_hotspots

    // TODO:
    //hotspot2.out.hotspots | footprints

    // Publish all our files
    Channel.empty().mix(
      bam_counts.out,
      count_adapters.out,
      dups.out.report,
      density.out.all,
      cutcounts.out,
      preseq.out,
      encode_cram.out.cram_with_index,
      spot_score.out.stats,
      insert_sizes.out,
      multimapping_density.out,
      normalize_density.out,
      motif_matrix.out,
    ).map { [it[0], params.outdir, it[1..-1].flatten()] } // Massage it into the right shape
    // Hotspots go in their own directory
    .mix(hotspot2.out.all.map { [it[0], "${params.outdir}/peaks", it[1..-1].flatten()] })
    | publish_many
}

process merge_bam {
  label "modules"
  // TODO: Optimzation for exactly one bam file

  input:
    // Inputs may be cram file but that's fine
    tuple val(meta), path('in*.bam')

  output:
    tuple val(meta), path('merged.bam')

  script:
    """
    samtools merge 'merged.bam' in*.bam
    """
}

// TODO: single end
process dups {
  label "modules"
  label 'high_mem'

  input:
    tuple val(meta), val(UMI), path(merged)

  output:
    tuple val(meta), file('marked.bam'), emit: bam
    tuple val(meta), file('MarkDuplicates.picard'), emit: report

  script:
    if (UMI) {
      cmd = "UmiAwareMarkDuplicatesWithMateCigar"
      extra = "UMI_TAG_NAME=XD"
    } else {
      cmd = "MarkDuplicatesWithMateCigar"
      extra = "MINIMUM_DISTANCE=300"
    }
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      "INPUT=${merged}" OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    picard "${cmd}" \
        INPUT=cigar.bam OUTPUT=marked.bam \
        $extra \
        METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'

    samtools index marked.bam
    """
}

process filter_bam {
  label "modules"

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), file("filtered.bam")

  script:
    flag = params.UMI ? 1536 : 512
    """
    samtools view -b -F "${flag}" marked.bam > filtered.bam
    """
}

process filter_nuclear {
  label "modules"
  input:
    tuple val(meta), path(nuclear_chroms), path(bam)

  output:
    tuple val(meta), file('nuclear.bam')

  script:
    """
    samtools index "${bam}"
    cat "${nuclear_chroms}" \
    | xargs samtools view -o nuclear.bam "${bam}"

    # debug
    samtools view -c nuclear.bam
    """
}

process macs2 {
  label "macs2"
  scratch false

  when:
    params.peakcaller == "macs2"

  input:
    tuple val(meta), file(bam)

  output:
    file 'NA_*'

  script:
    """
    macs2 callpeak \
      -t "$bam" \
      -f BAMPE \
      -g hs \
      -B \
      -q 0.01
    """
}

process hotspot2 {
  label "modules"

  container "fwip/hotspot2:latest"

  label 'high_mem'

  when:
    params.peakcaller == "hotspot2"

  input:
    tuple val(meta), val(hotspotid), path(mappable), path(chrom_sizes), path(centers), path(nuclear)


  output:
    tuple val(meta), file('peaks/nuclear*'), emit: all
    tuple val(meta), file('peaks/nuclear.hotspots.fdr0.05.starch'), emit: hotspots
    tuple val(meta), file('peaks/nuclear.peaks.fdr0.001.starch'), emit: peaks

  script:
    """
    export TMPDIR=\$(mktemp -d)
    hotspot2.sh -F 0.5 -p "varWidth_20_${hotspotid}" \
      -M "${mappable}" \
      -c "${chrom_sizes}" \
      -C "${centers}" \
      "${nuclear}" \
      'peaks'

    cd peaks

    # Rename peaks files to include FDR
    mv nuclear.peaks.narrowpeaks.starch nuclear.peaks.narrowpeaks.fdr0.05.starch
    mv nuclear.peaks.starch nuclear.peaks.fdr0.05.starch

    bash \$STAMPIPES/scripts/SPOT/info.sh \
      nuclear.hotspots.fdr0.05.starch hotspot2 nuclear.SPOT.txt \
      > nuclear.hotspot2.info

    # TODO: Move this to separate process
    hsmerge.sh -f 0.01 nuclear.allcalls.starch nuclear.hotspots.fdr0.01.starch
    hsmerge.sh -f 0.001 nuclear.allcalls.starch nuclear.hotspots.fdr0.001.starch

    density-peaks.bash \$TMPDIR "varWidth_20_${hotspotid}" nuclear.cutcounts.starch nuclear.hotspots.fdr0.01.starch ../"${chrom_sizes}" nuclear.density.starch nuclear.peaks.fdr0.01.starch \$(cat nuclear.cleavage.total)
    density-peaks.bash \$TMPDIR "varWidth_20_${hotspotid}" nuclear.cutcounts.starch nuclear.hotspots.fdr0.001.starch ../"${chrom_sizes}" nuclear.density.starch nuclear.peaks.fdr0.001.starch \$(cat nuclear.cleavage.total)

    rm -rf "\$TMPDIR"
    """
}

process spot_score {
  label "modules"

  input:
    tuple val(meta), val(genome_name), path(mappable), path(chrom_info), path(bam)

  output:
    tuple val(meta), file('r1.spot.out'), file('r1.hotspot.info'), emit: stats
    tuple val(meta), file('r1.spots.starch'), emit: spots

  script:
    """
    # random sample
    samtools view -h -F 12 -f 3 "$bam" \
      | awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
      | samtools view -uS -o __nuclear.bam
    bash \$STAMPIPES/scripts/bam/random_sample.sh __nuclear.bam subsample.bam 5000000
    samtools view -b -f 0x0040 subsample.bam > r1.bam

    # hotspot
    bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
      \$HOTSPOT_DIR \
      \$PWD \
      \$PWD/r1.bam \
      "${genome_name}" \
      "${params.readlength}" \
      DNaseI

    starch --header r1-both-passes/r1.hotspot.twopass.zscore.wig \
      > r1.spots.starch

    bash \$STAMPIPES/scripts/SPOT/info.sh \
      r1.spots.starch hotspot1 r1.spot.out \
      > r1.hotspot.info
    """
}

process bam_counts {
  label "modules"

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), file('tagcounts.txt')

  script:
    """
    python3 \$STAMPIPES/scripts/bwa/bamcounts.py \
      "$bam" \
      tagcounts.txt
    """
}

process count_adapters {
  label "modules"

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), file('adapter.counts.txt')

  script:
    """
    bash "\$STAMPIPES/scripts/bam/count_adapters.sh" "${bam}" \
    | sed 's/^/adapter\t/' \
    > adapter.counts.txt
    """
}

process preseq {
  label "modules"
  input:
    tuple val(meta), path(nuclear_bam)

  when:
    !params.UMI

  output:
    tuple val(meta), file('preseq.txt'), file('preseq_targets.txt'), file('dups.hist')

  script:
    """
    python3 \$STAMPIPES/scripts/bam/mark_dups.py -i "${nuclear_bam}" -o /dev/null --hist dups.hist
    preseq lc_extrap -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt \
    || preseq lc_extrap -defects -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt

    # write out preseq targets
    bash "\$STAMPIPES/scripts/utility/preseq_targets.sh" preseq.txt preseq_targets.txt
    """
}

process cutcounts {
  label "modules"
  label 'high_mem'


  input:
    tuple val(meta), path(fai), path(filtered_bam)

  output:
    tuple val(meta), file('fragments.starch'), file('cutcounts.starch'), file('cutcounts.bw'), file('cutcounts.bed.bgz'), file('cutcounts.bed.bgz.tbi')

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
    | starch - > cutcounts.starch

    # Bigwig
    "$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
      cutcounts.starch \
      cutcounts.bw \
      "${fai}"

    # tabix
    unstarch cutcounts.starch | bgzip > cutcounts.bed.bgz
    tabix -p bed cutcounts.bed.bgz
    """
}

process density {
  label "modules"

  label 'high_mem'

  input:
    tuple val(meta), path(chrom_bucket), path(fai), path(filtered_bam)

  output:
    tuple val(meta), file('density.starch'), file('density.bw'), file('density.bgz'), file('density.bgz.tbi'), emit: all
    tuple val(meta), file('density.starch'), emit: starch

  shell:
    window_size = 75
    bin_size = 20
    scale = 1_000_000
    '''
    mkfifo density.bed

    bam2bed -d \
    < "!{filtered_bam}" \
    | cut -f1-6 \
    | awk '{ if( $6=="+" ){ s=$2; e=$2+1 } else { s=$3-1; e=$3 } print $1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > density.bed \
    &

    unstarch "!{chrom_bucket}" \
    | bedmap --faster --echo --count --delim "\t" - density.bed \
    | awk -v "binI=!{bin_size}" -v "win=!{window_size}" \
          'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print $1 "\t" $2 + shiftFactor "\t" $3-shiftFactor "\tid\t" i $4}' \
    | starch - \
    > density.starch

    # Bigwig
    "$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
      density.starch \
      density.bw \
      "!{fai}" \
      "!{bin_size}"

    # Tabix
    unstarch density.starch | bgzip > density.bgz
    tabix -p bed density.bgz
    '''
}

process multimapping_density {
  label 'modules'
  label 'high_mem'

  input:
    tuple val(meta), path(chrom_bucket), path(fai), path(marked_bam)

  output:
    tuple val(meta), path("mm_density.starch"), path("mm_density.bw"), path('normalized.mm_density.starch'), path('normalized.mm_density.bw')

  shell:
    window_size = 75
    bin_size = 20
    scale = 1_000_000
    '''
    # Mark multi-mapping reads as QC-pass!
    samtools view -h "!{marked_bam}" |
    awk 'BEGIN{OFS="\t"} /XA:Z/ {$2 = and(or($2, 2), compl(512))} 1' |
    samtools view --threads 3 -F 512 -o filtered.bam
    samtools index filtered.bam


    # Generate density
    mkfifo density.bed

    bam2bed -d \
    < filtered.bam \
    | cut -f1-6 \
    | awk '{ if( $6=="+" ){ s=$2; e=$2+1 } else { s=$3-1; e=$3 } print $1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > density.bed \
    &

    unstarch "!{chrom_bucket}" \
    | bedmap --faster --echo --count --delim "\t" - density.bed \
    | awk -v "binI=!{bin_size}" -v "win=!{window_size}" \
          'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print $1 "\t" $2 + shiftFactor "\t" $3-shiftFactor "\tid\t" i $4}' \
    | starch - \
    > mm_density.starch

    # Bigwig
    "/home/solexa/stampipes/scripts/bwa/starch_to_bigwig.bash" \
      mm_density.starch \
      mm_density.bw \
      "!{fai}" \
      "!{bin_size}"

    # # Tabix
    # unstarch density.starch | bgzip > density.bgz
    # tabix -p bed density.bgz

    rm density.bed

    # Normalized density
    unstarch mm_density.starch \
      | awk -v allcounts=$(samtools view -c filtered.bam) \
            -v extranuclear_counts=$(samtools view -c "filtered.bam" chrM chrC) \
            -v scale=!{scale} \
            'BEGIN{ tagcount=allcounts-extranuclear_counts }
             { z=$5;
               n=(z/tagcount)*scale;
               print $1 "\t" $2 "\t" $3 "\t" $4 "\t" n }' \
      | starch - > normalized.mm_density.starch

    "$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
      normalized.mm_density.starch \
      normalized.mm_density.bw \
      "!{fai}" \
      "!{bin_size}"
    '''
}

process normalize_density {
  label "modules"

  input:
    tuple val(meta), path(fai), path(density), path(filtered_bam)

  output:
    tuple val(meta), file('normalized.density.starch'), file('normalized.density.bw'), file('normalized.density.bgz'), file('normalized.density.bgz.tbi')

  shell:
    bin_size = 20
    scale = 1_000_000
    '''
    samtools index "!{filtered_bam}"
    # Normalized density
    unstarch "!{density}" \
      | awk -v allcounts=$(samtools view -c !{filtered_bam}) \
            -v extranuclear_counts=$(samtools view -c "!{filtered_bam}" chrM chrC) \
            -v scale=!{scale} \
            'BEGIN{ tagcount=allcounts-extranuclear_counts }
             { z=$5;
               n=(z/tagcount)*scale;
               print $1 "\t" $2 "\t" $3 "\t" $4 "\t" n }' \
      | starch - > normalized.density.starch

    "$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
      normalized.density.starch \
      normalized.density.bw \
      "!{fai}" \
      "!{bin_size}"

    unstarch normalized.density.starch | bgzip > normalized.density.bgz
    tabix -p bed normalized.density.bgz
    '''
}

process insert_sizes {
  label "modules"

  input:
    tuple val(meta), path(nuclear_chroms), path(nuclear_bam)

  output:
    tuple val(meta), file('CollectInsertSizeMetrics.picard*')

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
  label "modules"

  input:
    tuple val(meta), path(fimo_transfac), path(fimo_names), path(hotspot_calls)
  //file fimo_transfac from file("${dataDir}/motifs/${genome_name}.fimo.starch")
  //file fimo_names from file("${dataDir}/motifs/${genome_name}.fimo.transfac.names.txt")

  output:
    tuple val(meta), file('hs_motifs*.txt')

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
  label "modules"

  input:
  tuple val(meta), path(transcript_starts), path(hotspot_calls)

  when:
    params.dofeatures

  output:
    tuple val(meta), file('prox_dist.info')

  script:
    thresholds = "0 1000 2500 5000 10000"
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

    for t in ${thresholds} ; do
      awk \
        -v t=\$t \
        '\$1 > t {sum+=1} END {print "percent-proximal-" t "bp " sum/NR}' \
        closest.clean.txt \
      >> prox_dist.info
    done
    """
}

/*
 * Metrics: Hotspot differential index comparison
 */
process differential_hotspots {
  label "modules"

  input:
    tuple val(meta), path(index), path(bam), path(peaks)

  output:
    tuple val(meta), file('differential_index_report.tsv')

  when:
    params.hotspot_index != "."

  shell:
    version = (new File(params.hotspot_index)).getAbsoluteFile().getParentFile().getName()
    diffName = "dhsindex_${version}_differential_peaks"
    diffPerName = "dhsindex_${version}_differential_peaks_percent"
    conName = "dhsindex_${version}_constitutive_peaks"
    conPerName = "dhsindex_${version}_constitutive_peaks_percent"

    '''
    set -e -o pipefail
    statOverlap=$(bedops -e 1 "!{peaks}" "!{index}" | wc -l)
    statNoOverlap=$(bedops -n 1 "!{peaks}" "!{index}" | wc -l)
    total=$(unstarch "!{peaks}" | wc -l)
    statPercOverlap=$(echo "scale=3; $statOverlap * 100.0/$total" | bc -q)
    statPercNoOverlap=$(echo "scale=3; $statNoOverlap * 100.0/$total" | bc -q)

    {
      echo -e "!{diffName}\t$statNoOverlap"
      echo -e "!{diffPerName}\t$statPercNoOverlap"
      echo -e "!{conName}\t$statOverlap"
      echo -e "!{conPerName}\t$statPercOverlap"
    } > differential_index_report.tsv
    '''
}


/*
 * Footprint calling
 */
// process learn_dispersion {
// 
//   label "footprints"
//   publishDir params.outdir
// 
//   memory = '8 GB'
//   cpus = 8
// 
//   when:
//     params.bias != ""
// 
//   input:
//     tuple val(meta), path(ref), path(bias), path(bam), path(spots)
// 
//   output:
//   set file('dm.json'), file(bam), file ("${bam}.bai") into dispersion
//   file 'dm.json' into to_plot
// 
//   script:
//   """
//   samtools index "$bam"
// 
//   # TODO: Use nuclear file
//   unstarch $spots \
//   | grep -v "_random" \
//   | grep -v "chrUn" \
//   | grep -v "chrM" \
//   > intervals.bed
// 
//   ftd-learn-dispersion-model \
//     --bm $bias \
//     --half-win-width 5 \
//     --processors 8 \
//     $bam \
//     $ref \
//     intervals.bed \
//   > dm.json
//   """.stripIndent()
// 
// }
// 
// process make_intervals {
// 
//   label "footprints"
//   input:
//   file starch from hotspot_calls
// 
//   output:
//   file 'chunk_*' into intervals mode flatten
// 
//   script:
//   """
//   unstarch "$starch" \
//   | grep -v "_random" \
//   | grep -v "chrUn" \
//   | grep -v "chrM" \
//   | split -l "$params.chunksize" -a 4 -d - chunk_
//   """.stripIndent()
// 
// }
// 
// process compute_deviation {
// 
//   label "footprints"
//   memory = '8 GB'
//   cpus = 4
// 
//   input:
//   set file(interval), file(dispersion), file(bam), file(bai) from intervals.combine(dispersion)
//   file(bias) from file(params.bias)
//   file(ref) from file("${params.genome}.fa")
// 
//   output:
//   file 'deviation.out' into deviations
// 
//   script:
//   """
//   ftd-compute-deviation \
//   --bm "$bias" \
//   --half-win-width 5 \
//   --smooth-half-win-width 50 \
//   --smooth-clip 0.01 \
//   --dm "$dispersion" \
//   --fdr-shuffle-n 50 \
//   --processors 4 \
//   "$bam" \
//   "$ref" \
//   "$interval" \
//   | sort --buffer-size=8G -k1,1 -k2,2n \
//   > deviation.out
//   """.stripIndent()
// }
// 
// process merge_deviation {
// 
//   label "footprints"
//   memory = "32 GB"
//   cpus = 1
// 
//   when:
//   params.bias != ""
// 
//   input:
//   file 'chunk_*' from deviations.collect()
// 
//   output:
//   file 'interval.all.bedgraph' into merged_interval
// 
//   script:
//   """
//   echo chunk_*
//   sort -k1,1 -k2,2n -S 32G -m chunk_* > interval.all.bedgraph
//   """.stripIndent()
// }
// 
// process working_tracks {
// 
//   label "footprints"
//   memory = '32 GB'
//   cpus = 1
// 
//   publishDir params.outdir
// 
//   input:
//   file merged_interval
// 
//   output:
//   file 'interval.all.bedgraph' into bedgraph
//   file 'interval.all.bedgraph.starch'
//   file 'interval.all.bedgraph.gz'
//   file 'interval.all.bedgraph.gz.tbi'
// 
//   script:
//   """
//   sort-bed "$merged_interval" | starch - > interval.all.bedgraph.starch
//   bgzip -c "$merged_interval" > interval.all.bedgraph.gz
//   tabix -0 -p bed interval.all.bedgraph.gz
//   """.stripIndent()
// 
// }
// 
// thresholds = Channel.from(0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
// 
// process compute_footprints {
// 
//   label "footprints"
//   memory = '8 GB'
//   cpus = 1
// 
//   publishDir params.outdir
// 
//   input:
//   set file(merged_interval), val(threshold) from merged_interval.combine(thresholds)
// 
//   output:
//   file "interval.all.fps.${threshold}.bed.gz"
//   file "interval.all.fps.${threshold}.bed.gz.tbi"
// 
//   script:
//   """
//   output=interval.all.fps.${threshold}.bed
//   cat "${merged_interval}" \
//   | awk -v OFS="\t" -v thresh="${threshold}" '\$8 <= thresh { print \$1, \$2-3, \$3+3; }' \
// 	| sort-bed --max-mem 8G - \
// 	| bedops -m - \
// 	| awk -v OFS="\t" -v thresh="${threshold}" '{ \$4="."; \$5=thresh; print; }' \
//   > \$output
// 
//   bgzip -c "\$output" > "\$output.gz"
//   tabix -0 -p bed "\$output.gz"
//   """.stripIndent()
// 
// }
// 
// process plot_footprints {
// 
//   label "footprints"
//   publishDir params.outdir
// 
//   input:
//   file model from to_plot
//   file plot from file("$baseDir/plot_footprints.py")
// 
//   output:
//   file "dispersion.*pdf"
// 
//   script:
//   """
//   "./$plot" "$model"
//   """
// }
