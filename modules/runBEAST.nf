outFolder = "${params.output}"

/*
 * Inject data into a specified BEAST xml template (fixed tree)
 */

process runBEAST {
    publishDir outFolder

    input:
        file beast_xml_file
        
    output:
        path("fixed_tree_epoch.log"), emit: logFile
        path("fixed_tree_epoch.trees"), emit: treesFile

    """
    beast ${beast_xml_file}
    """
}
