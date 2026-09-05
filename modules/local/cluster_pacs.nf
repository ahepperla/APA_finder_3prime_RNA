process CLUSTER_PACS {
    tag 'atlas-v1'
    label 'serial_high'

    publishDir "${params.outdir}/atlas", mode: 'copy', pattern: 'rejected_candidates.tsv.gz'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'pac_discovery.tsv'

    input:
    path evidence_tables
    path splice_continuations
    path normalized_samples
    path run_resolution
    path kernel
    path resolved_params

    output:
    path 'accepted_candidates.tsv.gz', emit: accepted
    path 'rejected_candidates.tsv.gz', emit: rejected
    path 'pac_discovery.tsv', emit: qc

    script:
    def knownArgument = params.known_pacs ? "--known-pacs '${params.known_pacs}'" : ''
    """
    pacusage cluster \
        --evidence ${evidence_tables.join(' ')} \
        --splice-continuations ${splice_continuations.join(' ')} \
        --samples '${normalized_samples}' \
        --resolution '${run_resolution}' \
        --kernel '${kernel}' \
        --params '${resolved_params}' \
        ${knownArgument} \
        --accepted accepted_candidates.tsv.gz \
        --rejected rejected_candidates.tsv.gz \
        --qc pac_discovery.tsv
    """
}
