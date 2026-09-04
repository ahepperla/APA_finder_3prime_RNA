include { CALIBRATE_LIBRARY_PROFILE } from '../../modules/local/calibrate_library_profile'
include { EXTRACT_3PRIME_EVIDENCE } from '../../modules/local/extract_3prime_evidence'
include { CLUSTER_PACS } from '../../modules/local/cluster_pacs'
include { ANNOTATE_PACS } from '../../modules/local/annotate_pacs'

workflow DISCOVERY {
    take:
    samples
    normalized_samples
    reference
    annotation
    resolved_params

    main:
    alignment_files = samples.map { meta, alignment, index, resolution, strand_qc, alignment_qc ->
        alignment
    }.collect()
    resolution_files = samples.map { meta, alignment, index, resolution, strand_qc, alignment_qc ->
        resolution
    }.collect()

    CALIBRATE_LIBRARY_PROFILE(
        normalized_samples,
        alignment_files,
        resolution_files,
        reference,
        annotation,
        resolved_params
    )
    run_resolution = CALIBRATE_LIBRARY_PROFILE.out.resolution.first()
    kernel = CALIBRATE_LIBRARY_PROFILE.out.kernel.first()

    EXTRACT_3PRIME_EVIDENCE(samples, reference, run_resolution, resolved_params)
    evidence_tables = EXTRACT_3PRIME_EVIDENCE.out.evidence
        .map { meta, tsv, parquet, plus, minus, qc -> tsv }
        .collect()
    CLUSTER_PACS(evidence_tables, run_resolution, kernel, resolved_params)
    ANNOTATE_PACS(
        CLUSTER_PACS.out.accepted,
        reference,
        annotation,
        run_resolution,
        resolved_params
    )

    emit:
    atlas = ANNOTATE_PACS.out.metadata
    atlas_bed = ANNOTATE_PACS.out.bed
    atlas_checksum = ANNOTATE_PACS.out.checksum
    motifs = ANNOTATE_PACS.out.motifs
    evidence = EXTRACT_3PRIME_EVIDENCE.out.evidence
    calibration = CALIBRATE_LIBRARY_PROFILE.out.calibration
    run_resolution = run_resolution
    kernel = kernel
    discovery_qc = CLUSTER_PACS.out.qc
    rejected = CLUSTER_PACS.out.rejected
    filtering_qc = EXTRACT_3PRIME_EVIDENCE.out.evidence.map {
        meta, tsv, parquet, plus, minus, qc -> qc
    }.collect()
}
