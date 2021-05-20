/*
 * Inject data into a specified BEAST xml template (fixed tree)
 */

process injectBEAST {
    input:
        file subtree_newick
        file aln_fa
        file dates_csv
        val clockrate
        
    output:
        fixed_tree.xml

    """
    beast_injector_fixed_tree_epoch.py --aln ${aln_fa} --dates ${dates_csv} --fixed_tree subtree_newick --clock_rate ${clockrate} --template ${template} > fixed_tree.xml
    """

}

