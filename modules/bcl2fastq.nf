nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

workflow test {
  main:
    def test_json = new JsonSlurper().parseText(file("/net/seq/data2/sequencers/220107_A01347_0081_BH7CW5DSX3/processing.json").text)
    def metadata = [
      "json": test_json,
      "flowcell_dir": file("/net/seq/data2/sequencers/220107_A01347_0081_BH7CW5DSX3"),
    ]
    Channel.of(metadata) | bcl2fastq
    bcl2fastq.out | view { it }
}

def group_bcl2fastq_output = { meta, fastq_files ->
  def out = []
  for (lib_data in meta.json.libraries) {
    def samplesheet_name = lib_data.samplesheet_name
    def lane = lib_data.lane
    def files = fastq_files.findAll { f ->
      def name = f.getName()
      name && name ==~ ~/${samplesheet_name}_S\d+_L00${lane}.*gz/
    }
    if (!files.isEmpty()) {
      out.add([lib_data, files.sort()])
    }
  }
  out
}

workflow bcl2fastq {

  take:
    // Format: ?
    metadata

  main:
    metadata.map { [it, it.json] } \
    | create_samplesheet

    // TODO take from metadata
    def num_lanes = 4

    create_samplesheet.out.combine( Channel.from(1..4) ) \
    | map { meta, samplesheet, lane -> [ meta, lane, samplesheet, meta.flowcell_dir ] } \
    | run_bcl2fastq

    out = run_bcl2fastq.out \
    | flatMap(group_bcl2fastq_output)

  emit:
    // Format: [meta, [r1.fastq.gz, r2.fastq.gz]]
    out
}

process create_samplesheet {
  input:
    tuple val(meta), val(json)

  output:
    tuple val(meta), file("Samplesheet.csv")

  script:
    def name = "Stamlab"
    def date = new java.text.SimpleDateFormat("MM/dd/yyyy").format(new Date())
    def header = """\
      [Header]
      Investigator Name,${name}
      Project Name,${name}
      Experiment Name,${name}
      Date,${date}
      Workflow,GenerateFASTQ
      [Settings]
      [Data]
      Lane,SampleID,SampleName,index,index2
      """.stripIndent()

    def cmd
    switch(json.flowcell.run_type) {
    case ~/.*Novaseq 6000.*/:
      cmd = """\
        jq -r '.libraries[] | select(.failed == false) | [(.lane|tostring), .samplesheet_name, .samplesheet_name, .barcode1.reverse_sequence, .barcode2.reverse_sequence,""] | join(",") ' \
        "processing.json" """.stripIndent()
      break
    default:
      throw new Exception("Unknown runtype ${json.flowcell.run_type}")
    }

    def processing_json = new JsonBuilder(json).toString()

    return """
cat  > processing.json <<EOF
${processing_json}
EOF
cat  > Samplesheet.csv <<EOF
${header}
EOF
${cmd} >> Samplesheet.csv
    """.stripIndent()
}

process run_bcl2fastq {

  cpus 20
  memory "32 GB"
  //module "bcl2fastq2"
  scratch false
  //container "genomicpariscentre/bcl2fastq2"
  container "altiusninstitute/stampipes-bcl2fastq2"

  input:
    tuple val(meta), val(lane), path("Samplesheet.csv"), path("flowcell_dir")

  output:
    tuple val(meta), file("output/*gz")

  script:
    // TODO
    def bcl_mask = meta.json.alignment_group.bases_mask
    def mismatches = 1
    """
    bcl2fastq \
      --runfolder-dir "\$(readlink flowcell_dir)" \
      --sample-sheet "Samplesheet.csv" \
      --use-bases-mask "${bcl_mask}" \
      --output-dir "output" \
      --barcode-mismatches "${mismatches}" \
      --tiles "s_${lane}_*"
    """
}
