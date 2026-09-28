process VALIDATE_INPUTS {
    tag 'inputs'
    label 'low'

    // RECORD_INPUT_CHECKSUMS and RECORD_SOFTWARE_VERSIONS complete and publish
    // input_checksums.tsv and software_versions.tsv.
    publishDir "${params.outdir}/manifest", mode: 'copy',
        pattern: '{resolved_params.yaml,normalized_samples.tsv,run_manifest.json}'
    publishDir "${params.outdir}/qc", mode: 'copy',
        pattern: '{input_validation.tsv,control_mapping.tsv}'

    input:
    path samplesheet
    val encoded_params
    path schema

    output:
    path 'resolved_params.yaml', emit: resolved_params
    path 'normalized_samples.tsv', emit: normalized_samples
    path 'software_versions.tsv', emit: software_versions
    path 'input_checksums.tsv', emit: input_checksums
    path 'run_manifest.json', emit: run_manifest
    path 'input_validation.tsv', emit: input_validation
    path 'control_mapping.tsv', emit: control_mapping

    script:
    """
    python -c 'import base64,sys; open("input_params.json","wb").write(base64.b64decode(sys.argv[1]))' '${encoded_params}'
    pacusage validate \
        --params input_params.json \
        --schema '${schema}' \
        --samples '${samplesheet}' \
        --output-dir .
    """
}

