include { QUANTIFY_PACS } from '../../modules/local/quantify_pacs'
include { MERGE_COUNTS } from '../../modules/local/merge_counts'

workflow QUANTIFICATION {
    take:
    evidence
    atlas
    run_resolution
    kernel
    resolved_params

    main:
    QUANTIFY_PACS(evidence, atlas, run_resolution, kernel, resolved_params)
    count_tables = QUANTIFY_PACS.out.counts.map { meta, table -> table }.collect()
    MERGE_COUNTS(count_tables, atlas, resolved_params)

    emit:
    wide = MERGE_COUNTS.out.wide
    long_counts = MERGE_COUNTS.out.long_counts
    gene_totals = MERGE_COUNTS.out.gene_totals
    pau = MERGE_COUNTS.out.pau
    testable = MERGE_COUNTS.out.testable
    filtering = MERGE_COUNTS.out.filtering
    quantification_qc = QUANTIFY_PACS.out.qc.collect()
}
