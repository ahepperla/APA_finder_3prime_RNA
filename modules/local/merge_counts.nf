process MERGE_COUNTS {
    tag 'all-samples'
    label 'medium'

    publishDir "${params.outdir}/counts", mode: 'copy',
        pattern: '{pac_counts.tsv.gz,pac_counts.long.parquet,gene_totals.tsv.gz,observed_pau.tsv.gz}'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'statistical_filtering.tsv'

    input:
    path count_tables
    path atlas
    path resolved_params

    output:
    path 'pac_counts.tsv.gz', emit: wide
    path 'pac_counts.long.parquet', emit: long_counts
    path 'gene_totals.tsv.gz', emit: gene_totals
    path 'observed_pau.tsv.gz', emit: pau
    path 'testable_pac_counts.tsv.gz', emit: testable
    path 'statistical_filtering.tsv', emit: filtering

    script:
    """
    pacusage merge-counts \
        --counts ${count_tables.join(' ')} \
        --atlas '${atlas}' \
        --wide pac_counts.tsv.gz \
        --long pac_counts.long.parquet \
        --gene-totals gene_totals.tsv.gz \
        --pau observed_pau.tsv.gz \
        --params '${resolved_params}' \
        --testable testable_pac_counts.tsv.gz \
        --filtering statistical_filtering.tsv
    """
}
