/*  
 * Perform ASR with gotree
 */
process gotreeASR {

    input:
        file newickFile
        file stateCsv

    output:
        path("asrtree.newick"), emit: asrtree

    """
    gotree acr --algo acctran -i ${newickFile} --states ${stateCsv}  -o asrtree.newick --random-resolve
    """
}

