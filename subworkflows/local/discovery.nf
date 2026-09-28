include { BUILD_CALIBRATION_REFERENCE } from '../../modules/local/build_calibration_reference'
include { SCAN_ALIGNMENT } from '../../modules/local/scan_alignment'
include { AGGREGATE_CALIBRATION } from '../../modules/local/aggregate_calibration'
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
    BUILD_CALIBRATION_REFERENCE(annotation)
    calibration_transcript_ends = BUILD_CALIBRATION_REFERENCE.out.transcript_ends.first()
    SCAN_ALIGNMENT(
        samples,
        reference,
        calibration_transcript_ends,
        resolved_params
    )
    calibration_summaries = SCAN_ALIGNMENT.out.calibration
        .map { meta, calibration -> calibration }
        .collect()
    AGGREGATE_CALIBRATION(
        normalized_samples,
        calibration_summaries,
        resolved_params
    )
    run_resolution = AGGREGATE_CALIBRATION.out.resolution.first()
    kernel = AGGREGATE_CALIBRATION.out.kernel.first()
    // A doubtful proximal-tag kernel is reported, not fatal: the run goes on.
    AGGREGATE_CALIBRATION.out.kernel_diagnostics
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.status == 'warning' }
        .subscribe { row ->
            log.warn "Calibration kernel: ${row.reason}. See qc/calibration_kernel_diagnostics.tsv."
        }

    EXTRACT_3PRIME_EVIDENCE(SCAN_ALIGNMENT.out.scan, run_resolution, resolved_params)
    evidence_tables = EXTRACT_3PRIME_EVIDENCE.out.evidence
        .map { meta, tsv, parquet, splice, plus, minus, qc -> parquet }
        .collect()
    splice_continuations = EXTRACT_3PRIME_EVIDENCE.out.evidence
        .map { meta, tsv, parquet, splice, plus, minus, qc -> splice }
        .collect()
    CLUSTER_PACS(
        evidence_tables,
        splice_continuations,
        normalized_samples,
        run_resolution,
        kernel,
        resolved_params
    )
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
    calibration = AGGREGATE_CALIBRATION.out.calibration
    kernel_diagnostics = AGGREGATE_CALIBRATION.out.kernel_diagnostics
    run_resolution = run_resolution
    kernel = kernel
    discovery_qc = CLUSTER_PACS.out.qc
    rejected = CLUSTER_PACS.out.rejected
    filtering_qc = EXTRACT_3PRIME_EVIDENCE.out.evidence.map {
        meta, tsv, parquet, splice, plus, minus, qc -> qc
    }.collect()
}
