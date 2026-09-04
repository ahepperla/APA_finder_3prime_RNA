process EXTRACT_3PRIME_EVIDENCE {
    tag { meta.sample_id }
    label 'high'

    publishDir "${params.outdir}/evidence", mode: 'copy',
        pattern: '*.{3prime_evidence.tsv.gz,3prime_evidence.parquet}'
    publishDir "${params.outdir}/tracks", mode: 'copy', pattern: '*.bedGraph.gz'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.fragment_filtering.tsv'

    input:
    tuple val(meta), path(alignment), path(index), path(sample_resolution),
        path(strandedness_qc), path(alignment_qc)
    path reference
    path run_resolution
    path resolved_params

    output:
    tuple val(meta),
        path("${meta.sample_id}.3prime_evidence.tsv.gz"),
        path("${meta.sample_id}.3prime_evidence.parquet"),
        path("${meta.sample_id}.plus.3prime_evidence.bedGraph.gz"),
        path("${meta.sample_id}.minus.3prime_evidence.bedGraph.gz"),
        path("${meta.sample_id}.fragment_filtering.tsv"),
        emit: evidence

    script:
    """
    pacusage extract-evidence \
        --sample-id '${meta.sample_id}' \
        --alignment '${alignment}' \
        --reference '${reference}' \
        --resolution '${run_resolution}' \
        --params '${resolved_params}' \
        --tsv '${meta.sample_id}.3prime_evidence.tsv.gz' \
        --parquet '${meta.sample_id}.3prime_evidence.parquet' \
        --plus-track '${meta.sample_id}.plus.3prime_evidence.bedGraph.gz' \
        --minus-track '${meta.sample_id}.minus.3prime_evidence.bedGraph.gz' \
        --qc '${meta.sample_id}.fragment_filtering.tsv' \
        --threads ${task.cpus}
    """
}

