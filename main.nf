#!/usr/bin/env nextflow

nextflow.enable.dsl=2

newickFile = file(params.newick)
alnFastaFile = file(params.aln)
dateCsvFile = file(params.csv)
outFolder = "${params.output}"

include { stateCsvFromNewick } from "./modules/stateCsvFromNewick.nf"
include { gotreeASR } from "./modules/gotreeASR.nf"
include { extractTransitionCsv } from "./modules/calSubtreeCsvFromASR.nf"
include { stratifiedSubsampleSequences } from "./modules/subsampleSequences.nf"
include { extractTreeSubset } from "./modules/extractTreeSubset.nf"
include { runTreeTime } from "./modules/runTreeTime.nf"
//include { injectBEAST } from "./modules/injectBEAST.nf"


process publish {

    publishDir outFolder

    input:
        file("sampleMetadataWithSubtree.csv")

    output:
        file("sampleMetadataWithSubtree.csv")

    """

    """

}

/*
 * Optional, slow version; perform ASR via MCMC on smaller subtree to assess CIs
 */
workflow asrMCMC {
    take:
        dateCsvFile
        alnFasta
        newick

    main:
        stratifiedSubsampleSequences(10000, dateCsvFile, alnFasta)
        extractTreeSubset(newick, stratifiedSubsampleSequences.out[0])
        runTreeTime(extractTreeSubset.out, stratifiedSubsampleSequences.out[1], stratifiedSubsampleSequences.out[2])
    
//    emit:


}

workflow {
    stateCsvFromNewick(newickFile)
    extractTransitionCsv(gotreeASR(newickFile, stateCsvFromNewick.out), dateCsvFile)
    publish(extractTransitionCsv.out[1])
    if (params.run_beast) {
        asrMCMC(dateCsvFile, alnFastaFile, newickFile)
    }
}
