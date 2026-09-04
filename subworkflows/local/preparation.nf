include { PREPARE_REFERENCE } from '../../modules/local/prepare_reference'
include { PREPARE_ALIGNMENT } from '../../modules/local/prepare_alignment'
include { INFER_STRANDEDNESS } from '../../modules/local/infer_strandedness'

workflow PREPARATION {
    take:
    normalized_samples
    resolved_params
    annotation

    main:
    PREPARE_REFERENCE(resolved_params)
    prepared_fasta = PREPARE_REFERENCE.out.fasta.first()

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

    PREPARE_ALIGNMENT(alignment_inputs, prepared_fasta)
    INFER_STRANDEDNESS(
        PREPARE_ALIGNMENT.out.prepared,
        prepared_fasta,
        annotation,
        resolved_params
    )

    emit:
    reference = prepared_fasta
    reference_index = PREPARE_REFERENCE.out.fai
    reference_qc = PREPARE_REFERENCE.out.metadata
    samples = INFER_STRANDEDNESS.out.resolved
    alignment_qc = PREPARE_ALIGNMENT.out.prepared.map {
        meta, alignment, index, qc -> qc
    }.collect()
    strandedness_qc = INFER_STRANDEDNESS.out.resolved.map {
        meta, alignment, index, resolution, qc, alignment_qc -> qc
    }.collect()
}
