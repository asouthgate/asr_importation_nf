/*
 * Date a tree with treetime
 */

process runTreeTime {
    input:
        file tree_newick
        file seq_fa
        file date_tsv

    output:
        file("timetree.newick")
    
    """
    sed 's/sequence_name/name/g' ${date_tsv} > names_dates.tsv
    treetime --tree ${tree_newick} --dates names_dates.tsv --aln ${seq_fa} --outdir .
    Rscript timetree.nexus timetree.newick
    """
}
