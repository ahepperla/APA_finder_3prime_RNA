process FINALIZE_USAGE_MODEL {
    tag "${family}"
    label 'serial_medium'

    input:
    tuple val(family), path(preliminary_directory), path(interval_files)

    output:
    tuple val(family), path('final-family-*'), emit: statistics

    script:
    def outputDirectory = "final-family-${task.index}"
    def intervalArguments = interval_files.join(',')
    """
    mkdir -p '${outputDirectory}'
    Rscript '${projectDir}/scripts/fit_usage_model.R' \
        --mode finalize \
        --preliminary-dir '${preliminary_directory}/preliminary' \
        --intervals '${intervalArguments}' \
        --output-dir '${outputDirectory}'
    """
}
