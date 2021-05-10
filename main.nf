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
        file("states.csv") into stateCsvChannel

    """
    cal_tiplab_state_mapping_from_tree.py ${newickFile} > states.csv
    """
}

/*  
 * Perform ASR with gotree
 */
process gotreeASR {

    publishDir outFolder

    input:
        file newickFile
        file stateCsv from stateCsvChannel

    output:
        file("asrtree.newick") into asrTreeChannel

    """
    gotree acr --algo acctran -i ${newickFile} --states ${stateCsv}  -o asrtree.newick --random-resolve
    """
}

/*  
 * Extract subtrees (taxon groups) that are putative independent
 *  imported clades.
 * Collect data together for each taxon to output final csv
 *  files with subtree designation and dates
 */
process extractTransitionCsv {

    publishDir outFolder

    input:
        file("asrtree.newick") from asrTreeChannel
        file dateCsvFile

    output:
        folder("root_output_folder"), file("imported_from.csv") into importedFromCsvChannel

    """
    extract_subtrees_state.py asrtree.newick -1 root_output_folder > imported_from.csv
    get_extracted_trees_csv.py root_output_folder/*[0-9] > clusters.csv
    merge_cluster_data_csv.py clusters.csv ${dateCsvFile} clusters_and_dates.csv
    """
}
