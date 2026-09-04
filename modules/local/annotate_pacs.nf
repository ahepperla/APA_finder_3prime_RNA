process ANNOTATE_PACS {
    tag 'atlas-v1'
    label 'high'

    publishDir "${params.outdir}/atlas", mode: 'copy', pattern: 'pacs.v1.*'
    publishDir "${params.outdir}/motifs", mode: 'copy', pattern: 'pac_motifs.tsv.gz'

    input:
    path candidates
    path reference
    path annotation
    path run_resolution
    path resolved_params

    output:
    path 'pacs.v1.metadata.tsv.gz', emit: metadata
    path 'pacs.v1.bed.gz', emit: bed
    path 'pacs.v1.sha256', emit: checksum
    path 'pac_motifs.tsv.gz', emit: motifs

    script:
    def knownArgument = params.known_pacs ? "--known-pacs '${params.known_pacs}'" : ''
    """
    pacusage annotate \
        --candidates '${candidates}' \
        --reference '${reference}' \
        --annotation '${annotation}' \
        --resolution '${run_resolution}' \
        --params '${resolved_params}' \
        ${knownArgument} \
        --metadata pacs.v1.metadata.tsv.gz \
        --bed pacs.v1.bed.gz \
        --motifs pac_motifs.tsv.gz \
        --checksum pacs.v1.sha256
    """
}
