process PREPARE_REFERENCE {
    tag 'reference'
    label 'low'
    // genome.fa links to the source FASTA; move keeps links as links, where
    // the default copy would duplicate the genome into the work directory.
    stageOutMode 'move'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'reference_preparation.tsv'
    publishDir "${params.outdir}/prepared_reference", mode: 'copy',
        pattern: 'genome.fa*', enabled: params.save_prepared_reference.toString() == 'true'

    input:
    path fasta, stageAs: 'source/*'

    // The FASTA and its index travel together, so no step rebuilds the index.
    output:
    tuple path('genome.fa'), path('genome.fa.fai'), emit: reference
    path 'reference_preparation.tsv', emit: metadata

    script:
    """
    pacusage prepare-reference \
        --fasta '${fasta}' \
        --output genome.fa \
        --metadata reference_preparation.tsv
    """
}
