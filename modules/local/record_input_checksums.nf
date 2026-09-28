process RECORD_INPUT_CHECKSUMS {
    tag 'inputs'
    label 'low'

    publishDir "${params.outdir}/manifest", mode: 'copy', pattern: 'input_checksums.tsv'

    // VALIDATE_INPUTS hashes the sample sheet, FASTA, and GTF. Each
    // PREPARE_ALIGNMENT task hashes its own alignment, in parallel.
    input:
    path validated_checksums, stageAs: 'validated/input_checksums.tsv'
    path normalized_samples
    path alignment_metadata

    output:
    path 'input_checksums.tsv', emit: checksums

    script:
    """
    pacusage record-input-checksums \
        --base '${validated_checksums}' \
        --samples '${normalized_samples}' \
        --alignment-metadata ${alignment_metadata} \
        --output input_checksums.tsv
    """
}
