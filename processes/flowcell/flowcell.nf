nextflow.enable.dsl = 2
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

params.help = false
params.metadata = ""
params.flowcell_dir = ""
params.out_dir = "output"

def helpMessage() {
  log.info"""
    Usage: nextflow run flowcell.nf \\
             --help
  """.stripIndent();
}
if (params.help) {
  helpMessage()
  exit(0)
}
dataDir = "${baseDir}/../../../data"


include { bcl2fastq } from "../../modules/bcl2fastq.nf"
include { publish_many } from "../../modules/utility.nf"


workflow {
  main:
    def processing_data = new JsonSlurper().parseText(file(params.metadata).text)
    def metadata = [
      "json": processing_data,
      "flowcell_dir": file(params.flowcell_dir),
    ]

    Channel.of(metadata) | bcl2fastq
    bcl2fastq.out
    | map { meta, files -> [meta, params.out_dir, files ]}
    | publish_many

}
