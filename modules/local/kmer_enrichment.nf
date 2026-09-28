process KMER_ENRICHMENT {
    tag 'exploratory-kmers'
    label 'serial_high'

    publishDir "${params.outdir}/motifs", mode: 'copy',
        saveAs: { value -> value.tokenize('/').last() }

    input:
    path statistics_files
    path atlas
    path resolved_params

    output:
    path 'kmer/*', emit: results

    script:
    """
    pacusage kmer-enrichment \
        --statistics . \
        --atlas '${atlas}' \
        --params '${resolved_params}' \
        --output-dir kmer
    """
}
