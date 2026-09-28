process PREPARE_ALIGNMENT {
    tag { meta.sample_id }
    label 'medium'
    // A sorted source is linked, not copied; move keeps links as links, where
    // the default copy would duplicate every BAM into the work directory.
    stageOutMode 'move'

    publishDir "${params.outdir}/qc", mode: 'copy',
        pattern: '*.alignment_preparation.tsv'
    publishDir "${params.outdir}/prepared_alignments", mode: 'copy',
        pattern: '*.{bam,bai,cram,crai}', enabled: params.save_prepared_alignments

    input:
    tuple val(meta), path(alignment)
    tuple path(reference), path(reference_index)

    output:
    tuple val(meta),
        path("${meta.sample_id}.${meta.extension}"),
        path("${meta.sample_id}.${meta.extension}.${meta.index_extension}"),
        path("${meta.sample_id}.alignment_preparation.tsv"),
        emit: prepared

    script:
    """
    pacusage prepare-alignment \
        --sample-id '${meta.sample_id}' \
        --alignment '${alignment}' \
        --reference '${reference}' \
        --output '${meta.sample_id}.${meta.extension}' \
        --metadata '${meta.sample_id}.alignment_preparation.tsv' \
        --threads ${task.cpus}
    """
}
