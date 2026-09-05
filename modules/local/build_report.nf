process BUILD_REPORT {
    tag 'report'
    label 'serial_high'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path artifacts

    output:
    path 'index.html', emit: html

    script:
    """
    pacusage report --results . --output index.html
    """
}
