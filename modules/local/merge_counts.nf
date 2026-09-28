process MERGE_COUNTS {
    tag 'all-samples'
    label 'serial_high'

    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path count_tables
    path atlas

    output:
    path 'pac_counts.tsv.gz', emit: wide
    path 'pac_counts.long.parquet', emit: long_counts
    path 'gene_totals.tsv.gz', emit: gene_totals
    path 'observed_pau.tsv.gz', emit: pau

    script:
    """
    pacusage merge-counts \
        --counts ${count_tables.join(' ')} \
        --atlas '${atlas}' \
        --wide pac_counts.tsv.gz \
        --long pac_counts.long.parquet \
        --gene-totals gene_totals.tsv.gz \
        --pau observed_pau.tsv.gz
    """
}
