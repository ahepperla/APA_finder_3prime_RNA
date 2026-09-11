process FIT_USAGE_MODEL {
    tag "${family}"
    label 'high'

    input:
    val family
    path normalized_samples
    path testable_counts
    path atlas
    path resolved_params
    path motif_scores
    path motif_sensitivity

    output:
    tuple val(family), path('family-*'), emit: statistics

    script:
    def outputDirectory = "family-${task.index}"
    """
    mkdir -p '${outputDirectory}'
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    Rscript '${projectDir}/scripts/fit_usage_model.R' \
        --family '${family}' \
        --samples '${normalized_samples}' \
        --counts '${testable_counts}' \
        --atlas '${atlas}' \
        --params '${resolved_params}' \
        --motif-scores '${motif_scores}' \
        --motif-sensitivity '${motif_sensitivity}' \
        --bootstrap-workers '${task.cpus}' \
        --output-dir '${outputDirectory}'
    """
}
