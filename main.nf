#!/usr/bin/env nextflow

newickFile = file(params.newick)
metadataCsvFile = file(params.csv)
outFolder = "${params.output}"

/* 
 * Convert tip labels to a state csv
 */
process tiplabs2statecsv {

}

/*  
 * Perform ASR with gotree
 */
process gotreeASR {

}

/*  
 * Extract subtrees that represent state transitions
 */
process extractTransitionSubtrees {

}

/*  
 * Merge subtree designations with date metadata
 */
process mergeSubtreeWithDates {

}


