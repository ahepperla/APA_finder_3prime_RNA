process FIT_USAGE_MODEL {
    tag 'all-comparisons'
    label 'serial_high'

    publishDir "${params.outdir}/statistics", mode: 'copy',
        pattern: '*.tsv.gz'
    publishDir "${params.outdir}/motifs", mode: 'copy',
        pattern: '*.preference*.tsv.gz'

    input:
    path normalized_samples
    path testable_counts
    path atlas
    path resolved_params
    path motif_scores
    path motif_sensitivity

    output:
    path '*.tsv.gz', emit: statistics

    script:
    """
    Rscript '${projectDir}/scripts/fit_usage_model.R' \
        --samples '${normalized_samples}' \
        --counts '${testable_counts}' \
        --atlas '${atlas}' \
        --params '${resolved_params}' \
        --motif-scores '${motif_scores}' \
        --motif-sensitivity '${motif_sensitivity}' \
        --output-dir .
    """
}
