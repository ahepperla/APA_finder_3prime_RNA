include { PREPARE_REFERENCE } from '../../modules/local/prepare_reference'
include { PREPARE_ALIGNMENT } from '../../modules/local/prepare_alignment'
include { INFER_STRANDEDNESS } from '../../modules/local/infer_strandedness'

workflow PREPARATION {
    take:
    normalized_samples
    resolved_params
    annotation

    main:
    PREPARE_REFERENCE(Channel.value(file(params.fasta, checkIfExists: true)))
    // The FASTA and its index travel together, so no step rebuilds the index.
    prepared_reference = PREPARE_REFERENCE.out.fasta
        .combine(PREPARE_REFERENCE.out.fai)
        .first()

    alignment_inputs = normalized_samples
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def extension = row.alignment.toLowerCase().endsWith('.cram') ? 'cram' : 'bam'
            def index_extension = extension == 'cram' ? 'crai' : 'bai'
            def meta = [
                sample_id: row.sample_id,
                condition: row.condition,
                control_condition: row.control_condition,
                layout: row.layout,
                strandedness: row.strandedness,
                library_profile: row.library_profile,
                evidence_source: row.evidence_source,
                extension: extension,
                index_extension: index_extension
            ]
            tuple(meta, file(row.alignment, checkIfExists: true))
        }

    PREPARE_ALIGNMENT(alignment_inputs, prepared_reference)
    INFER_STRANDEDNESS(
        PREPARE_ALIGNMENT.out.prepared,
        prepared_reference,
        annotation,
        resolved_params
    )
    prepared = PREPARE_ALIGNMENT.out.prepared.map { meta, alignment, index, qc ->
        tuple(meta.sample_id, meta, alignment, index, qc)
    }
    resolved = INFER_STRANDEDNESS.out.resolved.map { meta, resolution, qc ->
        tuple(meta.sample_id, resolution, qc)
    }
    samples = prepared
        .join(resolved, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, meta, alignment, index, alignment_qc, resolution, strandedness_qc ->
            tuple(meta, alignment, index, resolution, strandedness_qc, alignment_qc)
        }

    emit:
    reference = prepared_reference
    reference_qc = PREPARE_REFERENCE.out.metadata
    samples = samples
    alignment_qc = PREPARE_ALIGNMENT.out.prepared.map {
        meta, alignment, index, qc -> qc
    }.collect()
    strandedness_qc = INFER_STRANDEDNESS.out.resolved.map {
        meta, resolution, qc -> qc
    }.collect()
}
