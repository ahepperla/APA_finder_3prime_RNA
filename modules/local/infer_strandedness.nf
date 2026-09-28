process INFER_STRANDEDNESS {
    tag { meta.sample_id }
    label 'serial_medium'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.strandedness.tsv'

    input:
    tuple val(meta), path(alignment), path(index), path(alignment_qc)
    tuple path(reference), path(reference_index)
    path annotation
    path resolved_params

    // The alignment is not an output: passing it through made Nextflow copy
    // the whole file back from scratch. PREPARATION rejoins it by sample.
    output:
    tuple val(meta),
        path("${meta.sample_id}.resolution.json"),
        path("${meta.sample_id}.strandedness.tsv"),
        emit: resolved

    script:
    """
    pacusage infer-strandedness \
        --sample-id '${meta.sample_id}' \
        --alignment '${alignment}' \
        --alignment-metadata '${alignment_qc}' \
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
