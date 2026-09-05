process CALIBRATE_SAMPLE {
    tag { meta.sample_id }
    label 'high'

    input:
    tuple val(meta), path(alignment), path(index), path(sample_resolution),
        path(strandedness_qc), path(alignment_qc)
    path reference
    path transcript_ends
    path resolved_params

    output:
    tuple val(meta), path("${meta.sample_id}.calibration.json"), emit: calibration

    script:
    """
    pacusage calibrate-sample \
        --sample-id '${meta.sample_id}' \
        --alignment '${alignment}' \
        --resolution '${sample_resolution}' \
        --reference '${reference}' \
        --transcript-ends '${transcript_ends}' \
        --params '${resolved_params}' \
        --output '${meta.sample_id}.calibration.json' \
        --threads ${task.cpus}
    """
}
