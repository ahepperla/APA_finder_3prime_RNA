process CALIBRATE_LIBRARY_PROFILE {
    tag 'run'
    label 'high'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'library_calibration.tsv'
    publishDir "${params.outdir}/manifest", mode: 'copy',
        pattern: '{calibration_kernel.tsv,library_resolution.json}'

    input:
    path normalized_samples
    path alignments
    path resolutions
    path reference
    path annotation
    path resolved_params

    output:
    path 'library_calibration.tsv', emit: calibration
    path 'calibration_kernel.tsv', emit: kernel
    path 'library_resolution.json', emit: resolution

    script:
    """
    pacusage calibrate \
        --samples '${normalized_samples}' \
        --alignments ${alignments.join(' ')} \
        --resolutions ${resolutions.join(' ')} \
        --reference '${reference}' \
        --annotation '${annotation}' \
        --params '${resolved_params}' \
        --output library_calibration.tsv \
        --kernel calibration_kernel.tsv \
        --resolution library_resolution.json \
        --threads ${task.cpus}
    """
}

