process RECORD_SOFTWARE_VERSIONS {
    tag 'versions'
    label 'low'

    publishDir "${params.outdir}/manifest", mode: 'copy'

    // Loading the statistics and figure packages also stops a run that lacks
    // them in its first minutes. Kept apart from VALIDATE_INPUTS, so that
    // changing an R script reruns only the steps that use it on -resume.
    input:
    path python_versions, stageAs: 'python/software_versions.tsv'
    path statistics_script
    path figure_script

    output:
    path 'software_versions.tsv', emit: versions

    script:
    """
    cp '${python_versions}' software_versions.tsv
    Rscript '${statistics_script}' --mode versions --output software_versions.tsv
    Rscript '${figure_script}' --mode versions --output software_versions.tsv
    """
}
