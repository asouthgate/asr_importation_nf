#!/usr/bin/env nextflow

nextflow.enable.dsl=2

newickFile = file(params.newick)
dateCsvFile = file(params.csv)
outFolder = "${params.output}"

include { stateCsvFromNewick } from "./modules/stateCsvFromNewick.nf"
include { gotreeASR } from "./modules/gotreeASR.nf"
include { extractTransitionCsv } from "./modules/calSubtreeCsvFromASR.nf"

process publish {

    publishDir outFolder

    input:
        file("sampleMetadataWithSubtree.csv")

    output:
        file("sampleMetadataWithSubtree.csv")

    """

    """

}

workflow {
    stateCsvFromNewick(newickFile)
    extractTransitionCsv(gotreeASR(newickFile, stateCsvFromNewick.out), dateCsvFile)
    publish(extractTransitionCsv.out[1])
}
