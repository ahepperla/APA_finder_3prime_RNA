process INFER_STRANDEDNESS {
    tag { meta.sample_id }
    label 'serial_medium'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.strandedness.tsv'

    input:
    tuple val(meta), path(alignment), path(index), path(alignment_qc)
    path reference
    path annotation
    path resolved_params

    output:
    tuple val(meta), path(alignment), path(index),
        path("${meta.sample_id}.resolution.json"),
        path("${meta.sample_id}.strandedness.tsv"),
        path(alignment_qc),
        emit: resolved

    script:
    """
    pacusage infer-strandedness \
        --sample-id '${meta.sample_id}' \
        --alignment '${alignment}' \
        --reference '${reference}' \
        --annotation '${annotation}' \
        --profile '${meta.library_profile}' \
        --layout '${meta.layout}' \
        --strandedness '${meta.strandedness}' \
        --evidence-source '${meta.evidence_source}' \
        --endpoint-model '${params.endpoint_model}' \
        --minimum-informative ${params.strand_min_informative_fragments} \
        --maximum-sampled ${params.strand_max_sampled_fragments} \
        --decision-fraction ${params.strand_decision_fraction} \
        --random-seed ${params.random_seed} \
        --output '${meta.sample_id}.resolution.json' \
        --qc '${meta.sample_id}.strandedness.tsv'
    """
}
