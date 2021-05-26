maxCpu = params.max_cpu
outFolder = "${params.output}"

/*
 * Inject data into a specified BEAST xml template (fixed tree)
 */

process postprocessBEAST {

    publishDir outFolder

    input:
        file trees
        val root_date
        
    output:
        path("posterior_node_data.csv"), emit: posterior_node_csv
        path("binned_data.csv"), emit: binned_posterior_transition_csv
        
    """
    grep tree ${trees} > trees.noheader.newick
    get_posterior_node_data.py --trees trees.noheader.newick --nproc ${maxCpu} > posterior_node_data.csv
    get_binned_node_data.py --root_date ${root_date} --data posterior_node_data.csv > binned_data.csv
    """
}
