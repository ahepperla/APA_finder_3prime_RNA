process SCAN_ALIGNMENT {
    tag { meta.sample_id }
    label 'medium'

    input:
    tuple val(meta), path(alignment), path(index), path(sample_resolution),
        path(strandedness_qc), path(alignment_qc)
    tuple path(reference), path(reference_index)
    path transcript_ends
    path resolved_params

    // One read of the alignment gives calibration summaries for every
    // candidate evidence source and the evidence that EXTRACT writes.
    output:
    tuple val(meta), path("${meta.sample_id}.calibration.json"), emit: calibration
    tuple val(meta), path("${meta.sample_id}.scan.*"), emit: scan

    script:
    """
    pacusage scan-alignment \
        --sample-id '${meta.sample_id}' \
        --alignment '${alignment}' \
        --alignment-metadata '${alignment_qc}' \
        --resolution '${sample_resolution}' \
        --reference '${reference}' \
        --transcript-ends '${transcript_ends}' \
        --params '${resolved_params}' \
        --calibration '${meta.sample_id}.calibration.json' \
        --output-prefix '${meta.sample_id}.scan' \
        --threads ${task.cpus}
    """
}
