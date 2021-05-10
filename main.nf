#!/usr/bin/env nextflow

newickFile = file(params.newick)
dateCsvFile = file(params.csv)
outFolder = "${params.output}"

/* 
 * Parse states from newick tip labels
 * and output a csv mapping labels to states
 */
process stateCsvFromNewick {

    publishDir outFolder

    input:
        file newickFile

    output:
        file stateCsv into stateCsvChannel

    """
    cal_tiplab_state_mapping_from_tree.py ${newickFile} > stateCsv
    """
}

/*  
 * Perform ASR with gotree
 */
process gotreeASR {

    input:

    output:

    """

    """
}

/*  
 * Extract subtrees that represent state transitions
 */
process extractTransitionSubtrees {

    input:

    output:

    """

    """
}

/*  
 * Merge subtree designations with date metadata
 */
process mergeSubtreeWithDates {

    input:

    output:

    """

    """
}


