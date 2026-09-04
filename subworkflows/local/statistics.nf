include { MOTIF_SCORES } from '../../modules/local/motif_scores'
include { FIT_USAGE_MODEL } from '../../modules/local/fit_usage_model'
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
    FIT_USAGE_MODEL(
        normalized_samples,
        testable_counts,
        atlas,
        resolved_params,
        MOTIF_SCORES.out.primary,
        MOTIF_SCORES.out.sensitivity
    )
    KMER_ENRICHMENT(FIT_USAGE_MODEL.out.statistics, atlas)

    emit:
    results = FIT_USAGE_MODEL.out.statistics
    motif_scores = MOTIF_SCORES.out.primary
    motif_sensitivity = MOTIF_SCORES.out.sensitivity
    kmer_results = KMER_ENRICHMENT.out.results
}
