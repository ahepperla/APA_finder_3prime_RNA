process BOOTSTRAP_USAGE_INTERVALS {
    tag "${family}:${batch.baseName}"
    label 'statistics_bootstrap'
    maxForks (params.statistics_bootstrap_max_forks as int)

    input:
    tuple val(family), path(batch)
    path statistics_script

    output:
    tuple val(family), path('*.intervals.tsv.gz'), emit: intervals

    script:
    def outputName = "${batch.baseName}.intervals.tsv.gz"
    """
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    Rscript '${statistics_script}' \
        --mode bootstrap \
        --batch '${batch}' \
        --bootstrap-workers '${task.cpus}' \
        --output '${outputName}'
    """
}
