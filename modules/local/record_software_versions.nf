process RECORD_SOFTWARE_VERSIONS {
    tag 'versions'
    label 'low'

    publishDir "${params.outdir}/manifest", mode: 'copy'

    // Loading the statistics packages also stops a run that lacks them in its
    // first minutes. Kept apart from VALIDATE_INPUTS, so that changing the R
    // script reruns only the statistics on -resume.
    input:
    path python_versions, stageAs: 'python/software_versions.tsv'
    path statistics_script

    output:
    path 'software_versions.tsv', emit: versions

    script:
    """
    cp '${python_versions}' software_versions.tsv
    Rscript '${statistics_script}' --mode versions --output software_versions.tsv
    """
}
