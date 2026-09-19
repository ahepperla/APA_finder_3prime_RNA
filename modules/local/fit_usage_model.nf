process FIT_USAGE_MODEL {
    tag "${family}"
    label 'statistics_high'

    input:
    val family
    path normalized_samples
    path testable_counts
    path atlas
    path resolved_params
    path motif_scores
    path motif_sensitivity

    output:
    tuple val(family), path('family-*'), emit: preliminary
    tuple val(family), path('family-*/bootstrap-batches/*.rds'), emit: bootstrap_batches

    script:
    def outputDirectory = "family-${task.index}"
    """
    mkdir -p '${outputDirectory}'
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    Rscript '${projectDir}/scripts/fit_usage_model.R' \
        --mode fit \
        --family '${family}' \
        --samples '${normalized_samples}' \
        --counts '${testable_counts}' \
        --atlas '${atlas}' \
        --params '${resolved_params}' \
        --motif-scores '${motif_scores}' \
        --motif-sensitivity '${motif_sensitivity}' \
        --model-workers '${task.cpus}' \
        --bootstrap-batch-size '${params.statistics_bootstrap_batch_size}' \
        --output-dir '${outputDirectory}'
    """
}
