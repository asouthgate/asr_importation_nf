/*
 * Stratified subsampling of sequences by time
 */

process stratifiedSubsampleSequences {

    input:
        val n_sequences
        file metadata_csv

    output:
        file("metadata_sample.csv")
    
    """
    stratified_sample.py ${n_sequences} ${metadata_csv} metadata_sample.csv
    """
}
