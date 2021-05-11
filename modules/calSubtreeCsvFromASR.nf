/*  
 * Extract subtrees (taxon groups) that are putative independent
 *  imported clades.
 * Collect data together for each taxon to output final csv
 *  files with subtree designation and dates
 */
process extractTransitionCsv {

    input:
        file("asrtree.newick")
        file dateCsvFile

    output:
        tuple path("root_output_folder"), path("imported_from.csv")
        path("sample_metadata_with_clusters.csv")

    """
    mkdir root_output_folder
    extract_subtrees_state.py asrtree.newick -1 root_output_folder > imported_from.csv
    get_extracted_trees_csv.py root_output_folder/*[0-9] > clusters.csv
    merge_cluster_data_csv.py clusters.csv ${dateCsvFile} sample_metadata_with_clusters.csv
    """
}

