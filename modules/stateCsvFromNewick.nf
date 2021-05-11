/* 
 * Parse states from newick tip labels
 * and output a csv mapping labels to states
 */
process stateCsvFromNewick {

    input:
        file newickFile

    output:
        path("states.csv"), emit: state_csv

    """
    cal_tiplab_state_mapping_from_tree.py ${newickFile} > states.csv
    """
}

