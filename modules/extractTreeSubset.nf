/*
 * Extract a subset of a tree given a list of taxa names
 */

process extractTreeSubset {
    input:
        file newick
        file taxlabels
    output:
        file("subtree.newick")
    
    """
    extract_subtree_newick.R ${newick} ${taxlabels} subtree.newick
    """    
}
