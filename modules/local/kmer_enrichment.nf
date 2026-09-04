process KMER_ENRICHMENT {
    tag 'exploratory-kmers'
    label 'medium'

    publishDir "${params.outdir}/motifs", mode: 'copy',
        pattern: 'kmer/*.tsv.gz',
        saveAs: { value -> value.tokenize('/').last() }
    publishDir "${params.outdir}/motifs", mode: 'copy',
        pattern: 'kmer/kmer_enrichment_status.tsv',
        saveAs: { value -> value.tokenize('/').last() }

    input:
    path statistics_files
    path atlas

    output:
    path 'kmer/*', emit: results

    script:
    """
    mkdir -p kmer
    if [[ '${params.run_kmer_enrichment}' == 'true' ]]; then
        pacusage kmer-enrichment \
            --statistics . \
            --atlas '${atlas}' \
            --k ${params.motif_kmer_length} \
            --output-dir kmer
    else
        printf 'status\\nnot_requested\\n' > kmer/kmer_enrichment_status.tsv
    fi
    """
}
