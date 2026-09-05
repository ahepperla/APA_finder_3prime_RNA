process MERGE_USAGE_MODELS {
    tag 'all-families'
    label 'serial_medium'

    publishDir "${params.outdir}/statistics", mode: 'copy',
        pattern: 'merged/*.tsv.gz',
        saveAs: { value -> value.tokenize('/').last() }
    publishDir "${params.outdir}/motifs", mode: 'copy',
        pattern: 'merged/*.preference*.tsv.gz',
        saveAs: { value -> value.tokenize('/').last() }

    input:
    path family_directories

    output:
    path 'merged/*.tsv.gz', emit: statistics

    script:
    """
    pacusage merge-statistics \
        --inputs ${family_directories.join(' ')} \
        --output-dir merged
    """
}
