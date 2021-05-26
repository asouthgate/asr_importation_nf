outFolder = "${params.output}"

/*
 * Date a tree with treetime, returning a time tree and clock rate
 */

process runTreeTime {

    publishDir outFolder

    input:
        file tree_newick
        file seq_fa
        file date_tsv

    output:
        path "timetree.newick", emit: time_tree
        env CLOCKRATE, emit: clock_rate
        env ROOT_DATE, emit: root_date

    script:
    
    """
    echo ${date_tsv}
    sed 's/sequence_name/name/g' ${date_tsv} > names_dates.tsv
    treetime --tree ${tree_newick} --dates names_dates.tsv --aln ${seq_fa} --outdir . --clock-filter 3
    nexus2newick.R timetree.nexus timetree.newick
    CLOCKRATE=\$(python3 -c "with open('molecular_clock.txt') as f: print([l for l in f][1].split()[1])")
    ROOT_DATE=\$(sort -k2 dates.tsv | sed -n '2p' | cut -f2)
    """


}
