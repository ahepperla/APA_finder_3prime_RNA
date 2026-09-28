process MOTIF_SCORES {
    tag 'motif-usage'
    label 'serial_high'

    publishDir "${params.outdir}/motifs", mode: 'copy'

    input:
    path pau
    path atlas
    path resolved_params

    output:
    path 'motif_scores.tsv', emit: primary
    path 'motif_scores_known_rescue_sensitivity.tsv', emit: sensitivity

    script:
    """
    pacusage motif-scores \
        --pau '${pau}' \
        --atlas '${atlas}' \
        --params '${resolved_params}' \
        --output motif_scores.tsv
    pacusage motif-scores \
        --pau '${pau}' \
        --atlas '${atlas}' \
        --params '${resolved_params}' \
        --include-known-rescue \
        --output motif_scores_known_rescue_sensitivity.tsv
    """
}
