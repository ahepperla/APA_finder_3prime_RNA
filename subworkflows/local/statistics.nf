include { MOTIF_SCORES } from '../../modules/local/motif_scores'
include { FIT_USAGE_MODEL } from '../../modules/local/fit_usage_model'
include { BOOTSTRAP_USAGE_INTERVALS } from '../../modules/local/bootstrap_usage_intervals'
include { FINALIZE_USAGE_MODEL } from '../../modules/local/finalize_usage_model'
include { MERGE_USAGE_MODELS } from '../../modules/local/merge_usage_models'
include { KMER_ENRICHMENT } from '../../modules/local/kmer_enrichment'

workflow STATISTICS {
    take:
    normalized_samples
    testable_counts
    observed_pau
    atlas
    resolved_params

    main:
    MOTIF_SCORES(observed_pau, atlas)
    comparison_families = normalized_samples
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.condition != row.control_condition }
        .map { row -> row.control_condition }
        .unique()

    FIT_USAGE_MODEL(
        comparison_families,
        normalized_samples,
        testable_counts,
        atlas,
        resolved_params,
        MOTIF_SCORES.out.primary,
        MOTIF_SCORES.out.sensitivity
    )
    BOOTSTRAP_USAGE_INTERVALS(FIT_USAGE_MODEL.out.bootstrap_batches)
    bootstrap_intervals = BOOTSTRAP_USAGE_INTERVALS.out.intervals
        .groupTuple()
    finalization_inputs = FIT_USAGE_MODEL.out.preliminary
        .join(bootstrap_intervals)
    FINALIZE_USAGE_MODEL(finalization_inputs)
    family_directories = FINALIZE_USAGE_MODEL.out.statistics
        .map { family, directory -> directory }
        .collect()
    MERGE_USAGE_MODELS(family_directories)
    KMER_ENRICHMENT(MERGE_USAGE_MODELS.out.statistics, atlas)

    emit:
    results = MERGE_USAGE_MODELS.out.statistics
    motif_scores = MOTIF_SCORES.out.primary
    motif_sensitivity = MOTIF_SCORES.out.sensitivity
    kmer_results = KMER_ENRICHMENT.out.results
}
