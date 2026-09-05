process QUANTIFY_PACS {
    tag { meta.sample_id }
    label 'serial_medium'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.quantification.tsv'
    publishDir "${params.outdir}/counts/per_sample", mode: 'copy',
        pattern: '*.pac_counts.tsv.gz', enabled: params.save_intermediates

    input:
    tuple val(meta), path(evidence_tsv), path(evidence_parquet),
        path(splice_continuations), path(plus_track), path(minus_track), path(filter_qc)
    path atlas
    path run_resolution
    path kernel
    path resolved_params

    output:
    tuple val(meta), path("${meta.sample_id}.pac_counts.tsv.gz"), emit: counts
    path "${meta.sample_id}.quantification.tsv", emit: qc

    script:
    """
    pacusage quantify \
        --sample-id '${meta.sample_id}' \
        --evidence '${evidence_parquet}' \
        --atlas '${atlas}' \
        --resolution '${run_resolution}' \
        --kernel '${kernel}' \
        --params '${resolved_params}' \
        --counts '${meta.sample_id}.pac_counts.tsv.gz' \
        --qc '${meta.sample_id}.quantification.tsv'
    """
}
