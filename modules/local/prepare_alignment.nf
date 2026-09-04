process PREPARE_ALIGNMENT {
    tag { meta.sample_id }
    label 'medium'

    publishDir "${params.outdir}/qc", mode: 'copy',
        pattern: '*.alignment_preparation.tsv'
    publishDir "${params.outdir}/prepared_alignments", mode: 'copy',
        pattern: '*.{bam,bai,cram,crai}', enabled: params.save_prepared_alignments

    input:
    tuple val(meta), path(alignment)
    path reference

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

