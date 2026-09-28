process BUILD_REPORT {
    tag 'report'
    label 'serial_high'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path artifacts
    path resolved_params

    output:
    path 'index.html', emit: html

    script:
    """
    pacusage report --results . --params '${resolved_params}' --output index.html
    """
}
