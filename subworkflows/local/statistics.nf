include { MOTIF_SCORES } from '../../modules/local/motif_scores'
include { FIT_USAGE_MODEL } from '../../modules/local/fit_usage_model'
include { BOOTSTRAP_USAGE_INTERVALS } from '../../modules/local/bootstrap_usage_intervals'
include { FINALIZE_USAGE_MODEL } from '../../modules/local/finalize_usage_model'
include { MERGE_USAGE_MODELS } from '../../modules/local/merge_usage_models'
include { KMER_ENRICHMENT } from '../../modules/local/kmer_enrichment'

workflow STATISTICS {
    take:
    normalized_samples
    pac_counts
    observed_pau
    atlas
    resolved_params
    statistics_script

    main:
    MOTIF_SCORES(observed_pau, atlas, resolved_params)
    comparison_families = normalized_samples
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.condition != row.control_condition }
        .map { row -> row.control_condition }
        .unique()

    FIT_USAGE_MODEL(
        comparison_families,
        normalized_samples,
        pac_counts,
        atlas,
        resolved_params,
        MOTIF_SCORES.out.primary,
        MOTIF_SCORES.out.sensitivity,
        statistics_script
    )
    // One task per batch file; each family writes at least one batch.
    BOOTSTRAP_USAGE_INTERVALS(FIT_USAGE_MODEL.out.bootstrap_batches.transpose(), statistics_script)
    bootstrap_intervals = BOOTSTRAP_USAGE_INTERVALS.out.intervals
        .groupTuple()
    finalization_inputs = FIT_USAGE_MODEL.out.preliminary
        .join(bootstrap_intervals)
    FINALIZE_USAGE_MODEL(finalization_inputs, statistics_script)
    family_directories = FINALIZE_USAGE_MODEL.out.statistics
        .map { family, directory -> directory }
        .collect()
    MERGE_USAGE_MODELS(family_directories)
    KMER_ENRICHMENT(MERGE_USAGE_MODELS.out.statistics, atlas, resolved_params)

    emit:
    results = MERGE_USAGE_MODELS.out.statistics
    motif_scores = MOTIF_SCORES.out.primary
    motif_sensitivity = MOTIF_SCORES.out.sensitivity
    kmer_results = KMER_ENRICHMENT.out.results
}
