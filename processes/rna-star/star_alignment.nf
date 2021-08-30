nextflow.enable.dsl=2
params.outdir = "output"

params.r1 = null
params.r2 = null

params.adapter_p5 = null
params.adapter_p7 = null
params.readlength = null
params.umimethod = null

params.refdir = "/net/seq/data/genomes/human/GRCh38/noalts-sequins/"
starIndexDir = "${params.refdir}/STARgenome-gencode-v25"
params.star_threads = 8


include { star } from "./modules/star.nf" addParams(publish: false)
include { adapter_trim } from "../../modules/adapter_trimming.nf"
include { move_umt; takara_trim_umt } from "../../modules/umt.nf"
include { publish } from "../../modules/utility.nf"

workflow STAR_ALIGNMENT {
  
  main:
    ref_files = file("${starIndexDir}/*")
    adapter_trim( [params.r1, params.r2, params.adapter_p5, params.adapter_p7] )

    // Decide which UMI filtering to use, if any
    switch (params.umimethod.toLowerCase()) {
      case "takara-umt":
        takara_trim_umt(adapter_trim.out.fastq, params.readlength)
        star( takara_trim_umt.out.fastq, ref_files )
        move_umt(star.out.flatten().map {[it.getName(), it]})
        publish(move_umt.out)
        break

      case null:
      case "":
      case "none":
        star( adapter_trim.out.fastq, ref_files )
        publish(star.out.flatten().map {[it.getName(), it]})
        break

      default:
        throw new Exception("Can't understand umimethod parameter '${params.umimethod}'")
    }
}

workflow {
  STAR_ALIGNMENT()
}
