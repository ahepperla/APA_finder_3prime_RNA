process PREPARE_REFERENCE {
    tag 'reference'
    label 'medium'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'reference_preparation.tsv'
    publishDir "${params.outdir}/prepared_reference", mode: 'copy',
        pattern: 'genome.fa*', enabled: params.save_prepared_reference

    input:
    path resolved_params

    output:
    path 'genome.fa', emit: fasta
    path 'genome.fa.fai', emit: fai
    path 'reference_preparation.tsv', emit: metadata

    script:
    """
    FASTA=\$(python -c 'import yaml; print(yaml.safe_load(open("${resolved_params}"))["fasta"])')
    FAI=""
    if [[ -f "\${FASTA}.fai" ]]; then
        FAI="--fai \${FASTA}.fai"
    fi
    pacusage prepare-reference \
        --fasta "\${FASTA}" \
        \${FAI} \
        --output genome.fa \
        --metadata reference_preparation.tsv
    """
}

