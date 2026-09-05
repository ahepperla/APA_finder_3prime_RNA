process AGGREGATE_CALIBRATION {
    tag 'run'
    label 'low'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'library_calibration.tsv'
    publishDir "${params.outdir}/manifest", mode: 'copy',
        pattern: '{calibration_kernel.tsv,library_resolution.json}'

    input:
    path normalized_samples
    path calibrations
    path resolved_params

    output:
    path 'library_calibration.tsv', emit: calibration
    path 'calibration_kernel.tsv', emit: kernel
    path 'library_resolution.json', emit: resolution

    script:
    """
    pacusage aggregate-calibration \
        --samples '${normalized_samples}' \
        --calibrations ${calibrations.join(' ')} \
        --params '${resolved_params}' \
        --output library_calibration.tsv \
        --kernel calibration_kernel.tsv \
        --resolution library_resolution.json
    """
}
