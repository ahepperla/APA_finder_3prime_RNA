process PREPARE_REFERENCE {
    tag 'reference'
    label 'low'
    // genome.fa links to the source FASTA; move keeps links as links, where
    // the default copy would duplicate the genome into the work directory.
    stageOutMode 'move'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'reference_preparation.tsv'
    publishDir "${params.outdir}/prepared_reference", mode: 'copy',
        pattern: 'genome.fa*', enabled: params.save_prepared_reference

    input:
    path fasta, stageAs: 'source/*'

    output:
    path 'genome.fa', emit: fasta
    path 'genome.fa.fai', emit: fai
    path 'reference_preparation.tsv', emit: metadata

    script:
    """
    SOURCE=\$(python -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' '${fasta}')
    FAI=""
    if [[ -f "\${SOURCE}.fai" ]]; then
        FAI="--fai \${SOURCE}.fai"
    fi
    pacusage prepare-reference \
        --fasta '${fasta}' \
        \${FAI} \
        --output genome.fa \
        --metadata reference_preparation.tsv
    """
}
