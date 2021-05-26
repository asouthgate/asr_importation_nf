/*
 * Inject data into a specified BEAST xml template (fixed tree)
 */

process injectBEAST {
    input:
        file aln_fa
        file dates_csv
        file subtree_newick
        val clockrate
        
    output:
        path("fixed_tree.xml"), emit: beast_xml

    """
    beast_injector_fixed_tree_epoch.py --aln ${aln_fa} --dates ${dates_csv} --fixed_tree ${subtree_newick} --clock_rate ${clockrate} --xml_template template_fixed_tree_multiple_epochs.xml > fixed_tree.xml
    """

}

