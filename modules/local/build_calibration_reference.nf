process BUILD_CALIBRATION_REFERENCE {
    tag 'transcript-ends'
    label 'serial_medium'

    input:
    path annotation

    output:
    path 'calibration_transcript_ends.tsv', emit: transcript_ends

    script:
    """
    pacusage calibration-reference \
        --annotation '${annotation}' \
        --output calibration_transcript_ends.tsv
    """
}
