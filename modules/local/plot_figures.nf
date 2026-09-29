process PLOT_FIGURES {
    tag 'figures'
    label 'serial_medium'

    publishDir "${params.outdir}/figures", mode: 'copy',
        saveAs: { value -> value.tokenize('/').last() }

    input:
    path normalized_samples
    path statistics_files
    path atlas
    path observed_pau
    path resolved_params
    path figure_script

    output:
    path 'figures/*', emit: figures

    script:
    """
    Rscript '${figure_script}' \
        --mode figures \
        --samples '${normalized_samples}' \
        --statistics-dir . \
        --atlas '${atlas}' \
        --pau '${observed_pau}' \
        --params '${resolved_params}' \
        --output-dir figures
    """
}
