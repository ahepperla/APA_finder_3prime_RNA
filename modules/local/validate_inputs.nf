process VALIDATE_INPUTS {
    tag 'inputs'
    label 'low'

    // input_checksums.tsv is published by RECORD_INPUT_CHECKSUMS, which adds
    // the per-sample alignment hashes.
    publishDir "${params.outdir}/manifest", mode: 'copy',
        pattern: '{resolved_params.yaml,normalized_samples.tsv,software_versions.tsv,run_manifest.json}'
    publishDir "${params.outdir}/qc", mode: 'copy',
        pattern: '{input_validation.tsv,control_mapping.tsv}'

    input:
    path samplesheet
    val encoded_params

    output:
    path 'resolved_params.yaml', emit: resolved_params
    path 'normalized_samples.tsv', emit: normalized_samples
    path 'software_versions.tsv', emit: software_versions
    path 'input_checksums.tsv', emit: input_checksums
    path 'run_manifest.json', emit: run_manifest
    path 'input_validation.tsv', emit: input_validation
    path 'control_mapping.tsv', emit: control_mapping

    // The comment in the script records a change of output contract. Nextflow
    // hashes the script text, so it also stops -resume from reusing results
    // of older versions, whose input_checksums.tsv still had alignment rows.
    script:
    """
    python -c 'import base64,sys; open("input_params.json","wb").write(base64.b64decode(sys.argv[1]))' '${encoded_params}'
    # Alignments are hashed once each by PREPARE_ALIGNMENT, not here.
    pacusage validate \
        --params input_params.json \
        --samples '${samplesheet}' \
        --output-dir .
    """
}

