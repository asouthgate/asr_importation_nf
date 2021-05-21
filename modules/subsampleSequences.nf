/*
 * Stratified subsampling of taxa data by time
 */

process stratifiedSubsampleSequences {

    input:
        val n_sequences
        file metadata_csv
        file alignment_fasta

    output:
        file("sample_sequence_names.csv")
        file("sample.hfix.m.fa")
        file("sample_dates.tsv")
    
    """
    # Retrieve a stratified sample
    stratified_sample.py ${n_sequences} ${metadata_csv} metadata_sample.csv

    # Get the sequence names
    cut -f2 -d',' metadata_sample.csv > sample_sequence_names.csv

    # Fetch sequences
    extract_by_token_list.py ${alignment_fasta} sample_sequence_names.csv > sample.preqc.hfix.m.fa

    # Perform auto QC on the sequences
    all_qc.py sample.preqc.hfix.m.fa > sample.hfix.m.fa

    # Get metadata subsamples
    cut -f2,3 -d',' metadata_sample.csv > sample_dates.csv
    sed 's/,/\t/g' sample_dates.csv > sample_dates.tsv
    """
}
