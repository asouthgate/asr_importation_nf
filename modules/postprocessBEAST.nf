maxCpu = val(params.max_cpu)

/*
 * Inject data into a specified BEAST xml template (fixed tree)
 */

process runBEAST {
    input:
        file trees
        val root_date
        val max_date
        
    output:
        file("posterior_node_data.csv"), emit posterior_node_csv
        file("binned_data.csv"), emit binned_posterior_transition_csv
        
    """
    get_posterior_node_data.py --trees ${trees} --ncpu ${maxCpu} > posterior_node_data.csv
    get_binned_node_data.py --last_date ${max_date} --root_date ${root_date} --data posterior_node_data.csv > binned_data.csv
    """
}
